#!/usr/bin/env python
"""
Find evidence of co-residence between two sets of genomic features in a
metagenome assembly.

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
from seq_annot.seqio import open_io, write_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Beta"
__version__ = '0.1.1'


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

    def parse_gff(self, entry, tag='ID'):
        """Populate dictionary of cargo features by parsing GFF3Entry objects

        Args:
            entry (class): GFF3Entry object

            tag (str): GFF3 attribute to use as feature ID
        """
        ident = entry.attributes[tag]
        new_ident = "{}_{}".format(ident, len(self.cargo)+1)

        self.cargo[new_ident] = (entry.strand, entry.start, entry.end)

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


def output_dist(handle, contig, int1, int2):
    """Write contig length and relative distance distributions to file
    """
    if int1 and int2:
        idents2 = []
        for ident in int2:
            idents2 += contig.positional_ids(ident)

        for name in int1:
            # Find shortest distance between primary feature and any feature
            # in the secondary set
            dists = []
            idents1 = contig.positional_ids(name)
            for ident in idents1:
                dists.append(min([contig.distance(ident, k)[1] for k in \
                    idents2]))

            distance = min(dists)

            length = contig.length
            names2 = ";".join(sorted(int2))
            output = "{}\t{}\t{}\t{}\n".format(name, names2, length, distance)
            write_io(handle, output)

def increment_occurrence(co_res, contig, set1, set2, distribution=None):
    """Increment elements of the co-occurrence matrix when evidence of 
    co-redidence
    """
    features = set([i.rsplit("_", 1)[0] for i in contig.cargo.keys()])
    int1 = features.intersection(set1)
    int2 = features.intersection(set2)

    if distribution:
        output_dist(distribution, contig, int1, int2)

    row_pos = [set1.index(i) for i in int1] if int1 else [-1]
    col_pos = [set2.index(j) for j in int2] if int2 else [-1]
    for r in row_pos:
        for c in col_pos:
            co_res[r][c] += 1

    return(co_res)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('in_sets',
        metavar='in.csv',
        help="input feature sets file in CSV format. The sets file should "
             "contain exactly two columns of genomic feature names. Lines "
             "starting with '#' will be ignored. Use '-' to indicate that "
             "input should be taken from standard input (stdin)")
    parser.add_argument('in_gff',
        metavar='in.gff3',
        help="input feature annotation file in GFF3 format. Use '-' to indicate "
             "that input should be taken from standard input (stdin)")
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
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Argument sanity checks
    if (args.out_dist and not args.in_fasta):
        parser.error("error: argument -d/--out-dist must be supplied along "
            "with -f/--fasta")

    if args.in_sets == '-' and args.in_gff == '-':
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

    out_h = args.out 
    out_d = args.out_dist

    if args.in_gff == '-':
        in_gff = sys.stdin
    else:
        in_gff = args.in_gff

    if args.in_sets == '-':
        in_sets = sys.stdin
    else:
        in_sets = args.in_sets

    # Calculate and store contig lengths
    len_dist = {}
    if args.in_fasta:
        for record in fasta_iter(args.in_fasta):
            len_dist[record.id] = len(record.sequence)

    # Store elements of feature sets file
    set1 = []
    set2 = []
    with open_io(in_sets, mode='rb') as in_h:
        for nline, row in enumerate(in_h):
            row = row.decode('utf-8')

            if row.startswith('#'):  #skip comments
                continue

            row = row.strip().split('\t')

            try:
                el1, el2 = row
            except ValueError:
                raise FormatError("{}: line {}. Incorrect number of columns "
                    "provided. Please verify file format"\
                    .format(in_sets, nline))
            else:
                set1.append(el1)
                set2.append(el2)

    # Prepare array for storing occurrences of co-residence
    rows = set1 + ["None"]
    columns = set2 + ["None"]

    nrow = len(rows)
    ncol = len(columns)

    co_res = np.zeros((nrow, ncol))

    # Output header for distributions
    if out_d:
        write_io(out_d, "#ID\tColocated\tContigLength\tMinimumDistance\n")

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
        except KeyError:  #skip features without proper attribute tag
            no_attr += 1
            continue
         
        seqid = entry.seqid
        if seqid != giv.seqid:
            # Process previous genomic region
            co_res = increment_occurrence(co_res, giv, set1, set2, out_d) 

            # Reset for new genomic region interval
            giv = GenomicRegion()
            giv.seqid = seqid

            if out_d:  #add contig length to GenomicRegion object
                try:
                    giv.length = len_dist[seqid]
                except KeyError:
                    giv.length = "NA"

            ncontigs += 1
        
        giv.parse_gff(entry, attr_tag)  #add feature as cargo to contig

    # Process last contig
    co_res = increment_occurrence(co_res, giv, set1, set2, out_d)

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
