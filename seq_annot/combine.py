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
import os
from seq_annot.seqio import open_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.1.2'


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
    parser.add_argument('-p', '--precedence',
        choices=['order', 'none'],
        default='none',
        help="method to resolve overlapping features [default: none]. "
             "Options are 'order' and 'none'. If 'order', precedence will "
             "be determined by input file order")
    parser.add_argument('--preserve',
        dest='preserve_within',
        action='store_true',
        help="preserve features if they overlap within a GFF file")
    parser.add_argument('--stranded',
        dest='allow_diff_strand',
        action='store_true',
        help="preserve features if they overlap but only on different strands")
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
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables
    out_h = args.out.write
    out_d = args.discarded.write if args.discarded else do_nothing

    overlap_precedence = args.precedence

    gff_totals = 0
    o_totals = 0
    p_totals = 0

    uniques = {}  #store unique features
    chroms = {}  #store location of feature
    ind_totals = {}  #store individual GFF3 feature totals

    # Output unique features
    for position, gff in enumerate(args.gffs):
        gff_base = os.path.basename(gff)
        ind_totals[gff_base] = 0

        with open_io(gff) as gff_h:
            for feature in gff3_iter(gff_h):

                try:
                    seq_id = feature.seqid
                except AttributeError:
                    if feature.startswith('##'):  #headers
                        continue
                    elif feature.startswith('#'):
                        continue  #comments
                    else:
                        print("error: GFF3 file does not seem to be formatted "
                              "correctly", file=sys.stderr)
                        sys.exit(1)

                gff_totals += 1
                ind_totals[gff_base] += 1

                ident = feature.attributes["ID"]

                if "Parent" in feature.attributes:
                    parent_id = feature.attributes["Parent"]
                    try:

                        parent.add_child(feature)
                        continue
                    except AttributeError:
                        
                
                try:
                    chrom = uniques[seq_id]
                except KeyError:
                    p_totals += 1
                    uniques[seq_id] = [(position, feature)]
                    continue

                if overlap_precedence == "order":
                    is_unique = True  #assume feature is unique

                    start = feature.start
                    end = feature.end
                    strand = feature.strand

                    # Check if feature falls within the interval of another
                    for item in chrom:
                        chrom_feat = item[1]

                        # Allow overlap for features in the same GFF3 file
                        chrom_pos = item[0]
                        if position == chrom_pos and args.preserve_within:
                            continue

                        # Allow features to overlap on different strands
                        if args.allow_diff_strands and \
                            ((strand == '+' and chrom_feat.strand == '-') or \
                            (strand == '-' and chrom_feat.strand == '+')):
                            continue

                        iv_1 = set(range(chrom_feat.start, chrom_feat.end + 1))
                        iv_2 = set(range(start, end + 1))
                        if len(iv_1.intersection(iv_2)) > 0:
                            o_totals += 1
                            is_unique = False
                            break  #found overlap, no need to continue

                    if is_unique:
                        # Feature does not fall within the interval of another, so 
                        # add to uniques
                        p_totals += 1
                        uniques[seq_id].append((position, feature))
                    else:
                        out_d(feature.write().encode('utf-8'))
                else:
                    p_totals += 1
                    uniques[seq_id].append((position, feature))

    # Output combined GFF3
    header = "##gff-version 3\n"
    out_h(header.encode('utf-8'))

    for chrom in sorted(uniques):
        chrom_feature = 0
        for item in uniques[chrom]:
            chrom_feature += 1
            entry = item[1]

            # Update ID attribute to ensure uniqueness
            old_ident = entry.attributes['ID']
            new_ident = '{}_{!s}'.format(entry.seqid, chrom_feature)

            entry.attributes['ID'] = new_ident
            out_h(entry.write().encode('utf-8'))

            if entry.children:
                for child in entry.children:
                    chrom_feature += 1
                    entry.attributes['ID'] = '{}_{!s}'.format(child.seqid, chrom_feature)

                    # Update Parent attribute with new parent ID
                    child.attributes['Parent'] = new_ident
                    out_h(child.write().encode('utf-8'))

    # Calculate and print statistics
    print("Features processed:", file=sys.stderr)
    print("  - feature totals:\t{!s}".format(gff_totals), file=sys.stderr)
    for gff in ind_totals:
        print("    - features in {}:\t{!s}".format(gff, ind_totals[gff]), \
              file=sys.stderr)
    print("  - successfully merged:\t{!s}".format(p_totals), file=sys.stderr)
    print("  - discarded due to overlapping feature:\t{!s}".format(o_totals), \
          file=sys.stderr)
    print("", file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to combine {!s} features"\
          .format(total_time, gff_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
