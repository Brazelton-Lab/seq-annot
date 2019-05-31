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

import argparse
from bio_utils.iterators import fasta_iter, GFF3Reader
import logging
import numpy as np
from seq_annot.argparse import *
from seq_annot.seqio import *
import sys
import textwrap
from time import time

warn = logging.warning
logging.basicConfig(level=logging.ERROR)

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Beta"
__version__ = '0.4.3'

class FeatureNotFound(Exception):
    """A simple exception that is raised when an feature cannot be found on
    a genomic region
    """
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class GenomicRegion:
    """Class to store all feature annotations of a single contiguous genomic 
    region

    Attributes:
        seqid (str): ID of the genomic region (e.g. contiguous sequence, 
            scaffold, or genome)

        length (int): length of the genomic region in bp

        size (int): size of the genomic region in number of features predicted 
            along its length

        cargo (dict): dictionary of features contained within the genomic 
            region
    """

    def __init__(self):
        """Initialize variables to store information on the genomic region"""

        self.seqid = None
        self.length = None
        self.size = None
        self.cargo = {}

    def distance(self, f1, f2):
        """Calculate the physical distance in bp between two features 
        contained on the genomic region

        Args:
            f1 (str): name of the first feature contained within the genomic 
                region

            f2 (str): name of the second feature contained within the genomic 
                region
        
        Returns:
            int: physical distance in base-pairs between the two features
        """

        try:
            entry1 = self.cargo[f1]
            entry2 = self.cargo[f2]
        except KeyError as e:
            raise FeatureNotFound("feature {} not found on genomic region {}"\
                .format(e, self.seqid))

        if self.overlaps(f1, f2):  #overlapping features have zero distance
            return(0)

        iv1 = (entry1.start, entry1.end)
        iv2 = (entry2.start, entry2.end)

        # Calculate distance between features
        distance = min([abs(min(np.asarray(iv1) - i)) for i in iv2] + \
            [abs(min(np.asarray(iv2) - i)) for i in iv1])

        return(distance)

    def bounds(self, f1, f2, stranded: bool=False):
        """Determine if a genomic feature is bounded by another. A feature is 
        considered bounded if its interval (determined by the start and end 
        positions on the genomic region) is fully contained within the 
        interval of another feature

        Args:
            f1 (str): name of cargo feature possibly contained within the 
                second

            f2 (str): name of cargo feature possibly containing the first

            stranded (bool): a feature cannnot be bounded by another if on 
                different strands [default: False]

        Returns:
            bool: True if f1 is bounded by f2, else False
        """

        try:
            entry1 = self.cargo[f1]
            entry2 = self.cargo[f2]
        except KeyError as e:
            raise FeatureNotFound("feature {} not found on genomic region {}"\
                .format(e, self.seqid))

        # Determine if features reside on the same strand
        if stranded and ((entry1.strand == '+' and entry2.strand == '-') or \
            (entry1.strand == '-' and entry2.strand == '+')):
            return(False)  #feature is unbounded if on different strands

        iv1 = set(range(entry1.start, entry1.end + 1))
        iv2 = set(range(entry2.start, entry2.end + 1))

        return(iv1.issubset(iv2))

    def overlaps(self, f1, f2, stranded: bool=False):
        """Determine if the intervals of two features contained on a genomic 
        region overlap

        Args:
            f1 (str): name of feature in genomic region

            f2 (str): name of different feature in genomic region

            stranded (bool): don't allow features to be considered overlapping 
                if on different strands [default: False]

        Returns:
            bool: True if the intervals of f1 and f2 overlap, else False
        """

        try:
            entry1 = self.cargo[f1]
            entry2 = self.cargo[f2]
        except KeyError as e:
            raise FeatureNotFound("feature {} not found on genomic region {}"\
                .format(e, self.seqid))

        # Determine if features reside on the same strand
        if stranded and ((entry1.strand == '+' and entry2.strand == '-') or \
            (entry1.strand == '-' and entry2.strand == '+')):
            return(False)  #features dont overlap if on different strands

        iv1 = set(range(entry1.start, entry1.end + 1))
        iv2 = set(range(entry2.start, entry2.end + 1))

        return(len(iv1.intersection(iv2)) > 0)

    def intersect(self, names):
        """Find the intersection between the cargo features contained on the 
        genomic region and a set of feature names

        Args:
            names (list): list of feature names

        Returns:
            list: cargo feature IDs with corresponding elements in names
        """
        intersection = []
        # Populate list of cargo features found in set names
        for feature_id in self.cargo:
            feature = feature_id.split("_", 1)[0]
            if feature in names:
                intersection.append(feature_id)

        return(intersection)

    def add_cargo(self, entry, attr):
        """Populate dictionary of cargo features by parsing GFF3Entry objects

        Args:
            entry (class): GFF3Entry object

            attr (str): attribute to use as the feature identifier
        """
        try:
            feature = entry.attributes[attr]
        except KeyError:  #feature doesn't have proper attribute tag
            feature = "UNKNOWN"
            warn("record {} does not have attribute tag '{}'"\
                .format(entry.seqid, attr))

        # Create unique identifier for feature
        identifier = "{}_{}".format(feature, len(self.cargo) + 1)

        # Populate dictionary of cargo with new feature
        self.cargo[identifier] = entry


def output_dist(handle, giv, set1, set2, abunds={}, outall=False):
    """Write genomic region length and relative distance distributions to file

    Args:
        handle (file): file handle of output distributions file

        giv (class): GenomicRegion object containing co-occurring features

        set1 (list): list of elements in set1

        set2 (list): list of elements in set2

        abunds (dist): optional dictionary of abundances for a given genomic 
            feature

        outall (bool): output distributions for all genomic regions in GFF, not
            just for regions containing set elements [default: False]

    Return:
        exit code (int): returns 0 if able to write else 1
    """
    length = giv.length
    seqid = giv.seqid
    cargo = giv.cargo
    ncargo = len(cargo.keys())

    # Return 1 if unable to write to distributions file
    if not seqid:
        return(1)

    int1 = giv.intersect(set1)
    int2 = giv.intersect(set2)

    # Create dictionary of abundances. NAs if abunds not provided
    all_abunds = {}
    for element in int1 + int2:
        try:
            feature_id = cargo[element].attributes['ID']
        except KeyError:
            warn("Feature {} on genomic region {} does not contain the 'ID' "
                "attribute required for incorportaion of abundance"\
                .format(element, cargo[element].seqid))
            feature_id = 'NA'

        try:
            abund = abunds[feature_id]
        except KeyError:
            abund = 'NA'

        all_abunds[element] = abund

    # Handle case where completely overlapping sets and only one feature 
    # encountered with a corresponding element in the sets
    if int1 == int2 and len(int1) == 1:
        int2 = []

    if int1 and int2:  #contains elements from both sets
        used_ids = []
        for el1 in int1:
            name1, id1 = el1.split("_", 1)
            for el2 in int2:
                if el1 == el2:  #happens when sets intersection non-empty
                    continue

                name2, id2 = el2.split("_", 1)
                if id2 in used_ids:  #happens when sets intersection non-empty
                    continue

                # Find physical distance between elements
                distance = giv.distance(el1, el2)

                # Find abundance of elements
                abund1 = all_abunds[el1]
                abund2 = all_abunds[el2]

                bounded = 'TRUE' if giv.bounds(el1, el2) else 'FALSE'

                # Output information for given instance of co-occurrence
                output = "{}\t{}\t{!s}\t{!s}\t{}\t{!s}\t{!s}\t{!s}\t{}\n"\
                    .format(name1, name2, abund1, abund2, seqid, length, \
                    distance, ncargo, bounded)
                write_io(handle, output)

            used_ids.append(id1)

    elif int1 and not int2:  #contains elements only from first set
        for el1 in int1:
            name = el1.split("_", 1)[0]
            abund1 = all_abunds[el1]

            output = "{}\tNA\t{!s}\tNA\t{}\t{!s}\tNA\t{!s}\tNA\n".format(name,\
                abund1, seqid, length, ncargo)
            write_io(handle, output)

    elif int2 and not int1:  #contains elements only from secondary set
        for el2 in int2:
            name = el2.split("_", 1)[0]
            abund2 = all_abunds[el2]

            output = "NA\t{}\tNA\t{!s}\t{}\t{!s}\tNA\t{!s}\tNA\n".format(name,\
                abund2, seqid, length, ncargo)
            write_io(handle, output)

    else:
        if outall:
            output = "NA\tNA\tNA\tNA\t{}\t{!s}\tNA\t{!s}\tNA\n".format(seqid,\
                length, ncargo)
            write_io(handle, output)

    return(0)  #writes successful

def increment_occurrence(occur, giv, set1, set2, bounded=False, \
    symmetric=False):
    """Increment co-occurrence table

    Args:
        occur (Array): numpy Array containing incidences of co-occurrence

        giv (class): GenomicRegion object containing co-occurring features

        set1 (list): list of set1 elements

        set2 (list): list of set2 elements

        bounded (bool): in order to be considered co-occurring, a feature 
            found in the first set must be bounded by a feature within the 
            second set [default: False]

        symmetric (bool): generate symmetric matrix when identical elements 
            can be found in both sets [default: False]

    Returns:
        occur (Array): numpy Array incremented by the number instances of 
            co-occurrence
    """
    int1 = giv.intersect(set1) 
    int2 = giv.intersect(set2)

    if int1 and int2:  #contains elements from both sets
        used_ids = []
        for el1 in int1:
            name1, id1 = el1.split("_", 1)
            r = set1.index(name1)
            for el2 in int2:
                if el1 == el2:  #don't count same feature
                    continue

                name2, id2 = el2.split("_", 1)
                if id2 in used_ids:  #sets intersection non-empty
                    if not symmetric:
                        continue

                if bounded:  #el1 must be bounded by el2
                    # Determine if el2 bounds el1
                    if not giv.bounds(el1, el2):
                        continue

                c = set2.index(name2)

                occur[r][c] += 1

            used_ids.append(id1)

    return(occur)

def parse_sets(arg_in):
    """Parse input to sets arguments. If input is a file, read in elements of 
    the set line-by-line, otherwise assume set elements are separated by a 
    comma
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

def do_nothing(*args, **kwargs):
    return(1)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('in_gff',
        metavar='in.gff3',
        help="input feature annotations in GFF3 format. Use '-' to indicate "
             "that input should be taken from standard input (stdin)")
    parser.add_argument('-s', '--set',
        required=True,
        dest='in_set1',
        metavar='SET1',
        help="input primary feature sets as a file or comma-separated list. "
             "If input is a file, set elements should be line-separated. "
             "Use '-' to indicate that input should be taken from standard "
             "input (stdin)")
    parser.add_argument('-s2', '--set2',
        dest='in_set2',
        metavar='SET2',
        help="input optional supplementary feature sets as a file or "
            "comma-separated list. If input is a file, set elements should be "
            "line-separated. If unused, co-occurrences will be found between "
            "elements within the primary set")
    parser.add_argument('-f', '--fasta',
        dest='in_fasta',
        metavar='in.fa',
        help="input optional fasta file of contigs for calculation of contig "
             "lengths. When provided, contig lengths will be added to the "
             "output distributions file")
    parser.add_argument('-c', '--count',
        dest='in_count',
        metavar='in.csv',
        help="input optional tab-separated abundance file. The first column "
            "should consist of feature identifiers with matching values in the "
            "ID attribute of the GFF file. The second column should contain "
            "numerical values representing some metric of abundance for the "
            "corresponding feature. When provided, feature abundance will be "
            "added to the output distributions file")
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
        "just regions containing set elements [default: False]")
    parser.add_argument('--symmetric',
        action='store_true',
        help="output a symmetric co-occurrence table when the first and second "
            "sets have overlapping elements [default: False]")
    parser.add_argument('--bounded',
        action='store_true',
        help="increment co-occurrence table only if the interval of a feature "
            "found in the first set is fully bounded by the interval of a "
            "feature found in the second set [default: off]")
    parser.add_argument('--verbose',
        action='store_true',
        help="increase output verbosity")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Argument sanity checks
    if (args.in_fasta or args.counts) and not args.out_dist:
        parser.error("error: argument -d/--out-dist must be supplied whenever "
            "arguments -f/--fasta or -c/--count are used")

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
    sym = args.symmetric
    bounded = args.bounded

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

    if args.verbose:
        logging.basicConfig(level=logging.WARNING)

    # Calculate and store contig lengths
    len_dist = {}
    if args.in_fasta:
        with open_io(args.in_fasta, mode='rb') as fasta_h:
            for record in fasta_iter(fasta_h):
                len_dist[record.id] = len(record.sequence)

    # Store feature abundance values if provided
    abunds = {}
    if args.in_count:
        with open_io(args.in_count, mode='rb') as count_h:
            for nline, row in enumerate(count_h):

                row = row.decode('utf-8')

                if row.startswith('#'):  #skip comments
                    continue

                row = row.rstrip().split('\t')

                try:
                    name, value = row
                except ValueError:
                    raise FormatError("{}, line {!s}: Incorrect number of "
                        "columns provided. Please verify file format"\
                        .format(in_count, nline))

                if name not in abunds:
                    abunds[name] = value
                else:
                    raise FormatError("{}, line {!s}: duplicate feature IDs "
                        "encountered. IDs used in abundance file must be "
                        "unique".format(in_count, nline))

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
        header = "Set1\tSet2\tAbund1\tAbund2\tSeqID\tSeqLength\tMinDistance\t"\
            "NumFeatures\tBounded\n"
        write_io(out_d, header)

    # Iterate over GFF3 file
    nfeatures = 0
    no_attr = 0
    ncontigs = 1

    giv = GenomicRegion()
    gff_reader = GFF3Reader(open_io(in_gff, mode='rb'))
    for entry in gff_reader.iterate(parse_attr=True):
        nfeatures += 1

        seqid = entry.seqid
        old_seqid = giv.seqid
        if seqid != old_seqid:  #new contig encountered
            ncontigs += 1

            # Process previous genomic region
            co_res = increment_occurrence(co_res, giv, set1, set2, bounded, sym)

            # Output distributions
            ecode = distribution(out_d, giv, set1, set2, abunds, all_len)

            # Reset for new genomic region interval
            giv = GenomicRegion()
            giv.seqid = seqid

            # Add contig length to new GenomicRegion object
            try:
                giv.length = len_dist[seqid]
            except KeyError:
                giv.length = "NA"

        # Add feature as cargo to the genomic region
        giv.add_cargo(entry, attr_tag)

    # Process last contig
    co_res = increment_occurrence(co_res, giv, set1, set2, bounded, sym) 
    ecode = distribution(out_d, giv, set1, set2, abunds, all_len)

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
