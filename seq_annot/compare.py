#! /usr/bin/env python
"""
Compare annotated feature abundances across samples.

Usage:
    compare_features in.csv [in.csv ...]

Required input is one or more tabular files of feature abundances.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2).

Copyright:

    compare_features calculate feature abundance statistics
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

from arandomness.argparse import Open, ParseSeparator
import argparse
import csv
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
    parser.add_argument('csvs',
        metavar='in.csv [in.csv ...]',
        nargs='+',
        help="input one or more feature abundance files in CSV format")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output tabular feature by sample table [default: output to stdout]")
    parser.add_argument('-n', '--names',
        metavar='SAMPLE [,SAMPLE,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of sample names. The order should match "
             "the order of the input files")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()


    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('compare_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)


    # Track program run-time
    start_time = time()


    # Assign variables based on user input
    out_h = args.out.write

    sample_names = args.names if args.names else ["Sample{}".format(i) for i \
                   in list(range(len(args.csvs)))]


    # Store feature abundances in a dictionary
    totals = 0
    abundances = {}
    for position, csv_file in enumerate(args.csvs):
        with open_input(csv_file) as csv_h:
            for row in csv.reader(csv_h, delimiter='\t'):

                try:
                    name, abundance = row
                except:
                    raise FormatError("input file {} does not have the "
                                      "correct number of columns. Please "
                                      "verify that that the file is formatted "
                                      "correctly".format(csv_file))

                if name not in abundances:
                    abundances[name] = [0 for i in sample_names]

                abundances[name][position] = float(abundance)


    # Output header
    out_h("Feature\t{}\n".format('\t'.join(sample_names)))


    # Output feature abundances by sample
    for feature in sorted(abundances):
        out_h("{}\t{}\n".format(feature, '\t'.join([str(i) for i in \
              abundances[feature]])))


    # Output statistics
    feature_totals = len(abundances)
    sample_totals = len(sample_names)
    print("Total number of features:\t{!s}".format(feature_totals),\
          file=sys.stderr)
    print("Total samples merged:\t{!s}\n".format(sample_totals),\
          file=sys.stderr)


    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to merge {!s} samples\n"\
          .format(total_time, sample_totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
