#! /usr/bin/env python
"""
Combine GFF3 files, resolving any overlap conflicts through a hierarchy of 
precedence determined by input order.

Usage:
    combine_features [-o out.gff] in.gff [in.gff ...]

Required input is one or more GFF3 files of predicted features. Output is a 
GFF3 file containing combined features.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output 
to standard output (stdout).

Copyright:

    combine_features combine annotated features and resolve overlap conflicts
    Copyright (C) 2016  William Brazelton

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

from arandomness.argparse import Open
import argparse
from bio_utils.iterators import gff3_iter
from seq_annot.seqio import open_input
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.0.1'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gffs',
        metavar='in.gff [in.gff ...]',
        nargs='+',
        help="input one or more annotated feature files in GFF3 format")
    parser.add_argument('-o', '--out',
        metavar='out.gff',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output combined features in GFF3 format [default: output to "
             "stdout]")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('combine_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables
    out_h = args.out.write

    gff_totals = 0
    overlap_totals = 0
    passed_totals = 0

    uniques = {}  #to store unique features

    # Output unique features
    for gff in args.gffs:

        with open_input(gff) as gff_h:
            for entry in gff3_iter(gff_h):

                try:
                    seq_id = entry.seqid
                except AttributeError:
                    if entry.startswith('##'):  #output headers
                        out_h("{}\n".format(entry))
                        continue
                    else:
                        continue  #don't output comments

                start = entry.start
                end = entry.end
                strand = entry.strand
                
                gff_totals += 1

                try:
                    chrom = uniques[seq_id]
                except KeyError:
                    passed_totals += 1
                    entry.attributes['ID'] = 'f{!s}'.format(passed_totals)
                    uniques[seq_id] = [entry]  #first feature for given sequence
                    continue

                is_unique = True  #assume unique
                # Check if feature falls within the interval of another feature
                for chrom_feat in chrom:
                    # Allow features to overlap on different strands
                    if strand != chrom_feat.strand:
                        continue

                    interval = list(range(chrom_feat.start, chrom_feat.end + 1))
                    if start in int_range or end in interval:
                        overlap_totals += 1
                        is_unique = False
                        break  #found overlap, no need to continue

                if is_unique:
                    # Feature does not fall within the interval of another, so 
                    # add to uniques
                    passed_totals += 1
                    uniques[seq_id].append(entry)

    # Output combined GFF3
    for chrom in sorted(uniques):
        chrom_feature = 0
        for entry in chrom:
            chrom_feature += 1
            entry.attributes['ID'] = '{}_{!s}'.format(entry.seqid, chrom_feature)
            out_h(entry.write())

    # Calculate and print statistics
    print("Total features processed:\t{!s}".format(gff_totals), \
          file=sys.stderr)
    print("  - features remaining after merging:\t{!s}"\
          .format(passed_totals), file=sys.stderr)
    print("  - features overlapping a higher precedence feature:\t{!s}\n"\
          .format(overlap_totals), file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to combine {!s} features\n"\
          .format(total_time, gff_totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
