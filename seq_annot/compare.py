#! /usr/bin/env python
"""
Compare annotated features across samples.

Usage:
    compare_features in.csv [in.csv ...]

Required input is one or more tabular files of feature characteristics. Each 
input file should be composed of two or more columns, with the first column 
consisting of feature IDs. Subsequent columns should be numerical variables 
representing some form of information about a given feature, such as its 
abundance or prevalence within a sample. The default method for handling 
multiple variables is to sum the values (by row). Additional methods include 
taking the mean, median, and occurrence of each row or adding cell values.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2).

Copyright:

    compare_features compare feature statistics across samples
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
from collections import OrderedDict
import json
import numpy as np
import os
from seq_annot.seqio import open_io, write_io, FormatError
from statistics import mean, median
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.5.0'


def occurrence(values: list):
    """Find the number of items in a list with value greater than zero.

    Args:
        values (list): list of numerical values

    Retruns:
        int: number of items greater than zero
    """
    return len([i for i in values if i > 0])

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csvs',
        metavar='in.csv',
        nargs='+',
        help="input one or more feature characteristic files in CSV format")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output tabular feature x sample table [default: output to "
             "stdout]")
    parser.add_argument('-n', '--names',
        metavar='NAME [,NAME,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of names to be used as the table header "
             "[default: use input file basenames]. The order in which the "
             "names are given should correspond to the order of the input "
             "files if combining values using any method other than add")
    parser.add_argument('-i', '--ids',
        metavar='ID [,ID,...]',
        dest='features',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of feature IDs. Only those features "
             "in the list will be output [default: output all]")
    parser.add_argument('-m', '--method',
        metavar='METHOD',
        dest='method',
        default='sum',
        choices=["sum", "mean", "median", "occur", "add"],
        help="method for combining values, on a per feature basis, when input "
             "files contain more than two columns. Available options are add, "
             "sum, mean, median, and occur(ence) [default: sum]. If method is "
             "add, cells in the input tables will be added together. Adding "
             "cell values requires that the tables have equal dimensions. All "
             "other methods act on table rows")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Speedup trick
    list_type = type(list())

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('compare_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user input
    out_h = args.out
    infiles = args.csvs

    if args.method == "mean":
        method = mean
    elif args.method == "median":
        method = median
    elif args.method == "occur": 
        method = occurrence
    elif args.method == "add": 
        method = None
    else:
        method = sum

    colnames = args.names if args.names else [os.path.basename(i) for i \
                   in args.csvs]
    rownames = args.features

    # Combine feature abundances
    totals = 0
    abundances = OrderedDict()
    for position, infile in enumerate(infiles):
        with open_io(infile, mode='rb') as in_h:
            for nline, row in enumerate(in_h):

                row = row.decode('utf-8')

                if row.startswith('#'):  #skip comments
                    continue
                else:
                    totals += 1

                row = row.rstrip().split('\t')

                try:
                    name, values = row[0], row[1:]
                except ValueError:
                    raise FormatError("{}: line {}. Incorrect number of "
                        "columns provided. Please verify file format"\
                        .format(csv_file, nline))

                # Ignore features not found in list of features to include
                if rownames:
                    if name not in rownames:
                        continue

                try:
                    values = [float(j) for j in values]
                except ValueError:
                    raise FormatError("{}: line {}. Variables should only "
                        "contain numerical values".format(infile, nline))

                if method:
                    value = method(values)
                    try:
                        # Update sample value for existing feature
                        abundances[name][position] += value
                    except KeyError:
                        # Set default values for new feature and update current 
                        # sample value
                        abundances[name] = [float(0) for i in colnames]
                        abundances[name][position] += value

                else:
                    values = np.array(values)
                    try:
                        abundances[name] += values
                    except KeyError:
                        abundances[name] = values
                    except ValueError:
                        raise FormatError("{}: line {}. Input tables have "
                            "unequal dimensions".format(infile, nline))

    # Output header
    header = "#ID\t{}\n".format('\t'.join(colnames))
    write_io(out_h, header)

    # Output feature abundances by sample
    for rowname in abundances:
        row = abundances[rowname]
        row = '\t'.join(['{:g}'.format(i) for i in row])
        write_io(out_h, "{}\t{!s}\n".format(rowname, row))

    # Output statistics
    f_totals = len(abundances)
    s_totals = len(colnames)
    print("", file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - samples merged:\t{!s}".format(s_totals), file=sys.stderr)
    print("  - feature totals:\t{!s}".format(f_totals), file=sys.stderr)

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to merge {!s} features from {!s} samples"\
          .format(total_time, totals, s_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
