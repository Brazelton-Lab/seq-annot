#! /usr/bin/env python
"""
Combine abundance values of genomic features based on categories defined in an 
annotated GFF.

Required input is a GFF3 file of annotated protein-coding genes and a tabular 
file of feature abundances.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension to 
the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to 
standard output (stdout). Standard input (stdin) can be redirected to one of 
the positional argument if supplied with '-'.

Copyright:

    combine_abunds combines abundances using GFF to define categories
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

import argparse
from bio_utils.iterators import GFF3Reader
from collections import defaultdict
import logging
from seq_annot.argparse import *
from seq_annot.compare import occurrence
from seq_annot.seqio import *
from statistics import mean, median
import sys
import textwrap
from time import time

warn = logging.warning
logging.basicConfig(level=logging.ERROR)

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.0.7"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('in_gff',
        metavar='in.gff',
        help="input feature annotations in GFF3 format. Use '-' to indicate "
            "that input should be taken from standard input (stdin)")
    parser.add_argument('in_abund',
        metavar='in.abund',
        help="input tabular file of feature abundances. Use '-' to indicate "
            "that input should be taken from standard input (stdin)")
    parser.add_argument('-o', '--out',
        metavar='out.abund',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output tabular file of combined abundances [default: output "
            "to stdout]")
    parser.add_argument('-a', '--attr',
        metavar='TAG',
        dest='attr',
        help="attribute tag used to combine feature abundances [default: "
            "None]. If left off, abundances will be combined for a given "
            "genomic region")
    parser.add_argument('-i', '--id',
        metavar='TAG',
        dest='id',
        default='ID',
        help="attribute tag used to match GFF3 entries to feature IDs in the "
            "abundance file [default: ID]")
    parser.add_argument('-m', '--method',
        metavar='[sum|mean|median|occurrence]',
        dest='method',
        default='sum',
        choices=["sum", "mean", "median", "occurrence"],
        help="method used to combine feature abundance values [default: sum]")
    args = parser.parse_args()

    if args.in_abund == '-' and args.in_gff == '-':
        parser.error("error: standard input (stdin) can only be redirected to "
            "a single positional argument")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('combine_abunds', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Speedup tricks
    string_type = type(str())

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    if args.in_gff == '-':
        in_gff = sys.stdin
    else:
        in_gff = args.in_gff

    if args.in_abund == '-':
        in_abund = sys.stdin
    else:
        in_abund = args.in_abund

    out_h = args.out
    attr_tag = args.attr
    id_tag = args.id

    if args.method == "mean":
        method = mean
    elif args.method == "median":
        method = median
    elif args.method == "occur":
        method = occurrence
    else:
        method = sum

    # Iterate over abundance file, storing abundance values in dictionary
    totals = 0
    abunds = {}
    with open_io(in_abund, mode='rb') as in_h:
        for nline, row in enumerate(in_h):
            row = row.decode('utf-8')

            if row.startswith('#'):  #skip comments
                continue
            else:
                totals += 1

            row = row.rstrip().split('\t')

            try:
                name, value = row[0], row[1]
            except ValueError:
                infile = os.path.basename(args.in_abund)
                raise FormatError("{}, line {!s}: Incorrect number of "
                    "columns provided. Please verify file format"\
                    .format(infile, nline))

            try:
                value = float(value)
            except ValueError:
                infile = os.path.basename(args.in_abund)
                raise FormatError("{}, line {!s}: abundance value must be "
                    "either an integer or float".format(infile, nline))

            abunds[name] = value

    # Iterage over GFF, merging abundances of those features with the same 
    # attribute value
    nfeatures = 0
    cat_abunds = defaultdict(list)
    with open_io(in_gff, mode='rb') as gff_h:
        gff_reader = GFF3Reader(gff_h)
        for entry in gff_reader.iterate(parse_attr=True):
            nfeatures += 1

            try:
                feature_id = entry.attributes[id_tag]
            except KeyError:
                print("error: GFF entry {} does not have ID tag {} required "
                    "to link back to the abundance file"\
                    .format(entry.attributes['ID'], id_tag), file=sys.stderr)
                sys.exit(1)

            try:
                abund_value = abunds[feature_id]
            except KeyError:  #genomic feature has no associated abundance
                continue

            if attr_tag:
                try:
                    category = entry.attributes[attr_tag]
                except KeyError:
                    warn("GFF entry {} does not have attribute tag '{}'"\
                        .format(entry.attributes['ID'], attr_tag))
                    category = "UNKNOWN"
            else:
                category = entry.seqid  #abunds combined by genomic region

            cat_abunds[category].append(abund_value)

    # Output combined abundances
    for category in cat_abunds:
        final_value = method(cat_abunds[category])
        write_io(out_h, "{}\t{!s}\n".format(category, final_value))

    # Output statistics
    print("", file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - features combined:\t{!s}".format(totals), file=sys.stderr)
    print("  - output categories:\t{!s}".format(len(cat_abunds.keys())), \
        file=sys.stderr)
    print("  - GFF entries:\t{!s}".format(nfeatures), file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to merge abundances for {!s} features"\
          .format(total_time, totals), file=sys.stderr)
    print("", file=sys.stderr)

if __name__ == "__main__":
    main()
    sys.exit(0)
