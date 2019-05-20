#!/usr/bin/env python
"""
Find evidence of co-residence between two sets of genomic features in a
an assembly of mixed populations (such as a metagenome assembly).

Required inputs are a GFF3 file of annotated features and a CSV file 
containing exactly two columns. The first and second column should contain 
the content of the first and second set of genomic features, respectively. Set 
elements should have corresponding values within a single GFF attribute. 
Optional input is a FASTA file containing the contigs annotated in the GFF3 
file. Output is a co-occurrence table in CSV format.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2). Standard input (stdin) can be 
redirected to one of the positional arguments if supplied with '-'.

Copyright:

    colocate_features Finds instances of co-location between genomic features
    Copyright (C) 2018  William Brazelton

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from __future__ import print_function
from __future__ import division

from arandomness.argparse import Open,ParseSeparator
import argparse
from bio_utils.iterators import fasta_iter, GFF3Reader
import numpy as np
from seq_annot.seqio import open_io, write_io, FormatError
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Beta"
__version__ = '0.2.4'


class GenomicRegion:
    """Class to store all feature annotations of a single contiguous genomic 
    region

    Attributes:
        seqid (str): ID of sequence (e.g. contig)

        length (int): length of the genomic region in bp

        size (int): size of the genomic region in number of features predicted 
            to reside along its length

        cargo (dict): dictionary of features with their location along the 
            genomic region
    """

    def __init__(self):
        """Initialize variables to store information on the genomic region"""

        self.seqid = None
        self.length = None
        self.size = None
        self.cargo = {}

    def distance(self, f1, f2, stranded: bool=True):
        """Calculate the distance in bp between two features contained on the 
        genomic region

        Args:
            f1 (str): name of first feature

            f2 (str): name of second feature. Feature must already be

            stranded (bool): consider features encoded on different strands 
                [default: True]
        
        Returns:
            tuple: first element is the encoding strand (+/- if on separate 
                strands). Second element is an integer representing distance 
                between the features in bp
        """

        try:
            f1_strand, f1_start, f1_end = self.cargo[f1]
            f2_strand, f2_start, f2_end = self.cargo[f2]
        except KeyError:
            raise

        # Determine stranding of the two features
        if f1_strand == '+' and f2_strand == '+':
            strand = '+'
        elif f1_strand == '-' and f2_strand == '-':
            strand = '-'
        else:
            strand = '+/-'

        if not stranded and strand == '+/-':
            return((strand, "NA"))

        if self.overlap(f1, f2):  #overlapping features have zero distance
            return((strand, 0))

        iv_1 = (f1_start, f1_end)
        iv_2 = (f2_start, f2_end)

        # Calculate distance between features without regard to direction
        distance = min([abs(min(np.asarray(iv_1) - i)) for i in iv_2] + \
            [abs(min(np.asarray(iv_2) - i)) for i in iv_1])

        return((strand, distance))

    def overlap(self, f1, f2, stranded: bool=False):
        """Determine if two features in a genomic region overlap

        Args:
            f1 (str): name of feature in genomic region

            f2 (str): name of different feature in genomic region

            stranded (bool): don't allow features to be considered overlapping 
                if on different strands [default: False]

        Returns:
            bool: True if features overlap, else False
        """

        try:
            f1_strand, f1_start, f1_end = self.cargo[f1]
            f2_strand, f2_start, f2_end = self.cargo[f2]
        except KeyError:
            raise

        # Features can't be considered overlapping if on different strands
        if stranded and ((f1_strand == '+' and f2_strand == '-') or \
            (f1_strand == '-' and f2_strand == '+')):
            return(False)

        iv_1 = set(range(f1_start, f1_end + 1))
        iv_2 = set(range(f2_start, f2_end + 1))

        if len(iv_1.intersection(iv_2)) > 0:
            return(True)
        else:
            return(False)

    def bounded(self, f1, f2, stranded: bool=False):
        """Determine if a feature in a genomic region is bounded by 

        Args:
            f1 (str): name of feature in genomic region

            f2 (str): name of different feature in genomic region

            stranded (bool): a feature can only be bounded by others if on 
                different strands  [default: False]

        Returns:
            bool: True if features overlap, else False
        """

        pass

#        try:
#            f1_strand, f1_start, f1_end = self.cargo[f1]
#            f2_strand, f2_start, f2_end = self.cargo[f2]
#        except KeyError:
#            raise

        # Features can't be considered overlapping if on different strands
#        if stranded and ((f1_strand == '+' and f2_strand == '-') or \
#            (f1_strand == '-' and f2_strand == '+')):
#            return(False)


    def parse_gff(self, entry, feature):
        """Populate dictionary of cargo features by parsing GFF3Entry objects

        Args:
            entry (class): GFF3Entry object

            feature (str): feature identifier
        """
        pos_ident = "{}_{}".format(feature, len(self.cargo) + 1)

        self.cargo[pos_ident] = (entry.strand, entry.start, entry.end)

    def positional_ids(self, name):
        """Extract internal positional IDs

        Args:
            name (str): feature ID

        Returns:
            list: list of feature IDs with genomic position
        """
        ids = []
        for positional_id in self.cargo:
            if positional_id.rsplit("_", 1)[0] == name:
                ids.append(positional_id)

        return(ids)


def output_dist(handle, giv, set1, set2, outall=False):
    """Write genomic region length and relative distance distributions to file

    Args:
        handle (file): file handle of output distributions file

        giv (class): GenomicRegion object containing co-occurring features

        set1 (list): list of elements within set1

        set2 (list): list of elements within set2

        outall (bool): output distributions for all genomic regions in GFF, not
            just for regions containing set elements [default: False]

    Return:
        exit code (int): returns 0 if able to write else 1
    """
    length = giv.length
    seqid = giv.seqid
    cargo = giv.cargo.keys()
    ncargo = len(cargo)

    # Return 1 if unable to write to distributions file
    if not seqid:
        return(1)

    # Populate lists of features found in either or both sets
    idents1 = []
    idents2 = []
    for feature_id in cargo:
        feature = feature_id.split("_", 1)[0]
        if feature in set1:
            idents1.append(feature_id)
        if feature in set2:
            idents2.append(feature_id)

    if idents1 and idents2:  #contains elements from both sets
        used_ids = []
        for el1 in idents1:
            name1, id1 = el1.split("_", 1)
            for el2 in idents2:
                if el1 == el2:  #happens when sets intersection non-empty
                    continue

                name2, id2 = el2.split("_", 1)
                if id2 in used_ids:  #happens when sets intersection non-empty
                    continue

                # Find physical distance between elements
                strand, distance = giv.distance(el1, el2)

                # Output information for given instance of co-occurrence
                output = "{}\t{}\t{}\t{!s}\t{!s}\t{!s}\n".format(name1, \
                    name2, seqid, length, distance, ncargo)
                write_io(handle, output)

            used_ids.append(id1)

    elif idents1 and not idents2:  #contains elements only from first set
        for el1 in idents1:
            name = el1.split("_", 1)[0]
            output = "{}\tNA\t{}\t{!s}\tNA\t{!s}\n".format(name, seqid, \
                length, ncargo)
            write_io(handle, output)

    elif idents2 and not idents1:  #contains elements only from secondary set
        for el2 in idents2:
            name = el2.split("_", 1)[0]
            output = "NA\t{}\t{}\t{!s}\tNA\t{!s}\n".format(name, seqid, \
                length, ncargo)
            write_io(handle, output)

    else:
        if outall:
            output = "NA\tNA\t{}\t{!s}\tNA\t{!s}\n".format(seqid, length, \
                ncargo)
            write_io(handle, output)

    return(0)  #writes successful

def increment_occurrence(occur, giv, set1, set2):
    """Find instances of co-occurrence between set elements

    Args:
        occur (Array): numpy Array containing incidences of co-occurrence

        giv (class): GenomicRegion object containing co-occurring features

        set1 (list): list of elements within set1

        set2 (list): list of elements within set2

    Returns:
        occur (Array): numpy Array incremented by the number instances of 
            co-occurrence
    """
    
    # Convert feature positional IDs back to standard name
    features = [i.rsplit("_", 1)[0] for i in giv.cargo.keys()]

    # Find intersection between elements in set1 and features on the region
    int1 = set(features).intersection(set1)

    # Match intersection results to appropriate locations in co-occurrence table
    row_pos = []
    for i in int1:
        row_pos.append(set1.index(i))
        features.remove(i)

    # Find intersection between elements in set2 and remaining features
    int2 = set(features).intersection(set2)
    col_pos = [set2.index(j) for j in int2] if int2 else []

    # Increment co-occurrence table by row and column coordinates
    for r in row_pos:
        for c in col_pos:
            occur[r][c] += 1

    return(occur)

def parse_sets(arg_in):
    """Parse input to sets arguments. If input is a file, read in elements of 
    the set line-by-line, otherwise assume set elements are separated by a 
    comma.
    """
    in_set = []

    try:
        with open_io(arg_in, mode='rb') as in_h:
            for nline, row in enumerate(in_h):
                row = row.decode('utf-8')

                if row.startswith('#'):  #comments ignored
                    continue

                row = row.rstrip().split('\t')

                el = row[0]
                in_set.append(el)

    except FileNotFoundError:
        try:
            in_set = [i.lstrip() for i in arg_in.split(",")]
        except ValueError:
            raise ValueError("input provided to -s/--set or -s2/--set2 should "
                "either be a file or a comma-spearated string of set elements"\
                .format(arg_in))

    except ValueError:
        raise FormatError("{}: line {}. Incorrect number of columns provided. "
            "Please verify that the input file is correctly formatted"\
            .format(csv, nline))

    if not in_set:
        raise FormatError("something went wrong while attempting to parse "
            "input {}".format(arg_in))

    return(in_set)

def parse_set_entry():
    """
    """
    return(entry)

def do_nothing(*args, **kargs):
    return(1)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('in_gff',
        metavar='in.gff3',
        help="input feature annotation file in GFF3 format. Use '-' to indicate "
             "that input should be taken from standard input (stdin)")
    parser.add_argument('-s', '--set',
        required=True,
        dest='in_set1',
        metavar='in.csv',
        help="input primary feature sets file in CSV format. The file should "
             "contain a line-separated list of genomic feature names. Lines "
             "starting with '#' will be ignored. Use '-' to indicate that "
             "input should be taken from standard input (stdin)")
    parser.add_argument('-s2', '--set2',
        dest='in_set2',
        metavar='in.csv',
        help="input optional supplementary feature sets file in CSV format. "
             "The file should contain a line-separated list of genomic feature "
             "names. Lines starting starting with '#' will be ignored. If "
             "unused, co-occurrences will be found between elements within the "
             "primary set")
    parser.add_argument('-f', '--fasta',
        dest='in_fasta',
        metavar='in.fa',
        action=Open,
        mode='rb',
        help="input optional fasta file of contigs for calculation of contig "
             "lengths")
    parser.add_argument('-o', '--out',
        metavar='out.csv',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output co-occurrence table in CSV format [default: output to "
             "stdout]")
    parser.add_argument('-d', '--out-dist',
        dest="out_dist",
        metavar='out.csv',
        action=Open,
        mode='wb',
        help="output contig length and distance distributions in CSV format. "
             "Distributions are relative to the first set of features "
             "(primary set)")
    parser.add_argument('-a', '--attr',
        metavar='ATTRIBUTE',
        default="Alias",
        help="attribute tag for connecting GFF3 entries to elements within "
             "the feature sets file [default: Alias]")
    parser.add_argument('--all',
        action='store_true',
        help="output lengths for all genomic regions in the GFF3 file, not "
             "just regions containing set elements")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Argument sanity checks
    if (args.out_dist and not args.in_fasta):
        parser.error("error: argument -d/--out-dist must be supplied along "
            "with -f/--fasta")

    if args.in_set1 == '-' and args.in_gff == '-':
        parser.error("error: standard input (stdin) can only be redirected to "
            "a single positional argument")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('colocate_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    attr_tag = args.attr
    all_len = args.all

    out_h = args.out 
    out_d = args.out_dist

    if args.in_gff == '-':
        in_gff = sys.stdin
    else:
        in_gff = args.in_gff

    if args.in_set1 == '-':
        in_set1 = sys.stdin
    else:
        in_set1 = args.in_set1

    in_set2 = args.in_set2

    distribution = output_dist if out_d else do_nothing

    # Calculate and store contig lengths
    len_dist = {}
    if args.in_fasta:
        for record in fasta_iter(args.in_fasta):
            len_dist[record.id] = len(record.sequence)

    # Populate sets with elements from the CSV file
    set1 = parse_sets(in_set1)
    set2 = parse_sets(in_set2) if in_set2 else set1

    # Prepare data array for storing occurrences of co-residence
    rows = set1
    columns = set2

    nrow = len(rows)
    ncol = len(columns)

    co_res = np.zeros((nrow, ncol))

    # Output header for distributions file
    if out_d:
        header = "Set1\tSet2\tSeqID\tSeqLength\tMinDistance\tNumFeatures\n"
        write_io(out_d, header)

    # Iterate over GFF3 file
    nfeatures = 0
    no_attr = 0
    ncontigs = 1

    giv = GenomicRegion()
    gff_reader = GFF3Reader(open_io(in_gff, mode='rb'))
    for entry in gff_reader.iterate(parse_attr=True):
        nfeatures += 1

        try:
            feature = entry.attributes[attr_tag]
        except KeyError:  #feature doesn't have proper attribute tag
            no_attr += 1
            if all_len:
                feature = "UNKNOWN"
            else:
                continue
         
        seqid = entry.seqid
        old_seqid = giv.seqid
        if seqid != old_seqid:  #new contig encountered

            # Process previous genomic region
            co_res = increment_occurrence(co_res, giv, set1, set2) 

            # Output distributions
            ecode = distribution(out_d, giv, set1, set2, all_len)

            # Reset for new genomic region interval
            giv = GenomicRegion()
            giv.seqid = seqid

            # Add contig length to GenomicRegion object
            if out_d:
                try:
                    giv.length = len_dist[seqid]
                except KeyError:
                    giv.length = "NA"

            ncontigs += 1

        # Add feature as cargo to the genomic region
        giv.parse_gff(entry, feature)

    # Process last contig
    co_res = increment_occurrence(co_res, giv, set1, set2)

    ecode = distribution(out_d, giv, set1, set2, all_len)

    # Output co-occurrence table
    write_io(out_h, "ID\t{}\n".format("\t".join(columns)))
    for row_num, row_name in enumerate(rows):
        output = "{}\t{!s}\n".format(row_name, \
            "\t".join([str(i) for i in co_res[row_num]]))
        write_io(out_h, output)
    
    # Output statistics
    print("", file=sys.stderr)
    print("Contigs processed:", file=sys.stderr)
    print("  - contig totals:\t{!s}".format(ncontigs), file=sys.stderr)
    print("  - feature totals:\t{!s}".format(nfeatures), file=sys.stderr)
    print("  - features without attribute:\t{!s}".format(no_attr), file=sys.stderr)
    print("  - total instances of co-occurrence:\t{!s}".format(co_res.sum()), \
        file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to search for co-resident features on {!s} "
          "contigs".format(total_time, ncontigs), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
   main()
   sys.exit(0)
