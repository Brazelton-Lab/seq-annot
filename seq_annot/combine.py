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
from bio_utils.iterators import GFF3Reader
import os
from seq_annot.seqio import open_io, write_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.2.0'


def do_nothing(args):
    pass


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
        mode='wb',
        default=sys.stdout,
        help="output combined features in GFF3 format [default: output to "
             "stdout]")
    parser.add_argument('--precedence',
        dest='precedence',
        action='store_true',
        help="resolve feature overlap using input file order to determine "
             "feature precedence")
    parser.add_argument('--preserve',
        dest='preserve_within',
        action='store_true',
        help="preserve features if they overlap within a GFF file")
    parser.add_argument('--stranded',
        dest='stranded',
        action='store_true',
        help="preserve features if they overlap, but only on different strands")
    parser.add_argument('-d', '--discarded',
        metavar='out.gff',
        action=Open,
        mode='wb',
        help="output features in GFF3 format discarded due to overlapping "
             "intervals")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('combine_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables
    out_h = args.out
    out_d = args.discarded.write if args.discarded else do_nothing

    no_overlaps = args.precedence
    allow_diff_strands = args.stranded

    gff_totals = 0
    o_totals = 0
    p_totals = 0
    ind_totals = {}  #store individual GFF3 feature totals

    uniques = {}  #store unique features
    file_order = {}
    # Output unique features
    for file_number, gff in enumerate(args.gffs):
        gff_base = os.path.basename(gff)
        ind_totals[gff_base] = 0
        file_order[file_number] = gff_base

        with open_io(gff) as gff_h:
            gff_reader = GFF3Reader(gff_h)

            for feature in gff_reader.iterate():
                try:
                    chrom_id = feature.seqid
                except AttributeError:
                    if feature.startswith('##'):  #headers
                        continue
                    elif feature.startswith('#'):
                        continue  #comments
                    else:
                        line_number = gff_reader.current_line
                        print("error: {}, line {!s}: unable to parse GFF3 "
                              "file. The file may be formatted incorrectly"\
                              .format(gff_base, line_number), file=sys.stderr)
                        sys.exit(1)

                gff_totals += 1
                ind_totals[gff_base] += 1

                try:
                    chrom = uniques[chrom_id]
                except KeyError:  # First feature encountered at given pos
                    p_totals += 1
                    uniques[chrom_id] = [(file_number, feature)]
                    continue

                if no_overlaps:
                    is_unique = True  #assume feature is unique

                    # Check if feature falls within the interval of another
                    for item in chrom:
                        chrom_feat = item[1]

                        # Allow overlap for features in the same GFF3 file
                        chrom_pos = item[0]
                        if file_number == chrom_pos and args.preserve_within:
                            continue

                        if feature.overlap(chrom_feat):
                            o_totals += 1
                            is_unique = False
                            break  #found overlap, no need to continue

                    if is_unique:
                        # Feature does not fall within the interval of another, so 
                        # add to uniques
                        p_totals += 1
                        uniques[chrom_id].append((file_number, feature))
                    else:
                        out_d(feature.write().encode('utf-8'))
                else:
                    p_totals += 1
                    uniques[chrom_id].append((file_number, feature))

    # Output combined GFF3
    header = "##gff-version 3\n"
    write_io(out_h, header)

    id_lookup = {}
    for chrom in sorted(uniques):
        chrom_feature = 0
        for item in uniques[chrom]:
            chrom_feature += 1
            entry = item[1]

            write_io(out_h, entry.write())

    # Calculate and print statistics
    print("", file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - feature totals:\t{!s}".format(gff_totals), file=sys.stderr)
    print("  - features merged:\t{!s}".format(p_totals), file=sys.stderr)
    if args.precedence:
        print("  - overlaps discarded:\t{!s}"\
              .format(o_totals), file=sys.stderr)
    print("GFFs processed:", file=sys.stderr)
    print("  - files combined:\t{!s}".format(len(args.gffs)), file=sys.stderr)
    for gff in ind_totals:
        print("    - features in {}:\t{!s}".format(gff, ind_totals[gff]), \
              file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to combine {!s} features from {!s} files"\
          .format(total_time, gff_totals, len(args.gffs)), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
