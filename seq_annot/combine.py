#! /usr/bin/env python
"""
Combine GFF3 files, resolving any overlap conflicts either through a hierarchy 
of precedence determined by input order or the score column of the GFF files.

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

import argparse
from bio_utils.iterators import GFF3Reader
import os
from seq_annot.argparse import *
from seq_annot.seqio import *
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.3.3'


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
    parser.add_argument('-d', '--discards',
        metavar='out.gff',
        action=Open,
        mode='wb',
        help="output features in GFF3 format discarded due to overlapping "
             "intervals")
    parser.add_argument('-c', '--conflict',
        metavar='[order|score]',
        choices=["score", "order"],
        help="method used to resolve conflicting feature annotations [default: "
            "score]")
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
    out_d = args.discards.write if args.discards else do_nothing
    res_method = args.conflict

    gff_totals = 0
    o_totals = 0
    p_totals = 0
    ind_totals = {}  #store individual GFF3 feature totals

    regions = {}  #for storing features residing on genomic regions
    file_order = {}
    # Output unique features
    for file_number, gff in enumerate(args.gffs):
        gff_base = os.path.basename(gff)
        ind_totals[gff_base] = 0
        file_order[file_number] = gff_base

        with open_io(gff) as gff_h:
            gff_reader = GFF3Reader(gff_h)

            for entry in gff_reader.iterate():
                try:
                    chrom_id = entry.seqid
                except AttributeError:
                    if entry.startswith('##'):  #headers
                        continue
                    elif entry.startswith('#'):
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
                    ident = entry.attributes['ID']
                except KeyError:
                    line_number = gff_reader.current_line
                    raise FormatError("{}, line {!s}: ID is a required GFF "
                        "attribute".format(gff_base, line_number))

                try:
                    chrom = regions[chrom_id]
                except KeyError:  # First time encountering genomic region
                    p_totals += 1
                    regions[chrom_id] = GenomicRegion()
                    regions[chrom_id].seqid = chrom_id
                    regions[chrom_id].cargo[ident] = entry
                    continue
                else:
                    cargo = regions[chrom_id].cargo

                if ident in cargo:
                    o_totals += 1
                    if res_method == "order":  #feature already annotated
                        continue
                    else:
                        # Resolve conflicting annotation using score column
                        score = entry.score
                        prev_score = cargo[ident].score
                        if score <= prev_score:  #existing entry is better
                            continue
                        else:  #existing entry worse than current entry
                            cargo[ident] = entry
                else:
                    cargo[ident] = entry

    # Output combined GFF3
    header = "##gff-version 3\n"
    write_io(out_h, header)

    for chrom_id in sorted(regions):
        chrom = regions[chrom_id]
        for feature_id in chrom.order():
            cargo = chrom.cargo[feature_id]

            write_io(out_h, cargo.write())

    # Calculate and print statistics
    print("", file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - genomic regions:\t{!s}".format(p_totals), file=sys.stderr)
    print("  - entries encountered:\t{!s}".format(gff_totals), file=sys.stderr)
    print("  - conflicts discarded:\t{!s}".format(o_totals), file=sys.stderr)
    print("GFFs processed:", file=sys.stderr)
    print("  - files combined:\t{!s}".format(len(args.gffs)), file=sys.stderr)
    for gff in ind_totals:
        print("    - entries in {}:\t{!s}".format(gff, ind_totals[gff]), \
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
