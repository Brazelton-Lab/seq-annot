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
import json
import os
from seq_annot import merge_dbs
from seq_annot.seqio import open_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.2.1'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csvs',
        metavar='in.csv',
        nargs='+',
        help="input one or more feature abundance files in CSV format")
    parser.add_argument('-m', '--mapping',
        metavar='in.json',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input one or more relational databases in JSON format containing "
             "supplementary information about the features that should be "
             "added to the table")
    parser.add_argument('-f', '--fields',
        metavar='FIELD [,FIELD,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database "
            "that should be added to the table")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output tabular feature by sample table [default: output to stdout]")
    parser.add_argument('-n', '--names',
        metavar='DATASET [,DATASET,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of dataset names to be used as the header "
             "[default: will use file names as dataset names]. The order "
             "should match the order of the input files")
    parser.add_argument('-i', '--ids',
        metavar='FEATURE [,FEATURE,...]',
        dest='features',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of feature IDs. Only those features "
             "with a match in the list will be output [default: output all]")
    derep_group = parser.add_mutually_exclusive_group()
    derep_group.add_argument('-r', '--dup-file',
        metavar='in.tsv',
        dest='dup_file',
        action=Open,
        mode='rt',
        help="input tab-separated file containing database entries to combine. "
             "The input relational databases should be provided as columns and "
             "entries to merge between databases should be provided as rows. "
             "Lines starting with # will be ignored.")
    derep_group.add_argument('-e', '--dup-field',
        metavar='FIELD',
        dest='dup_field',
        help="database field by which entries should be combined. The content "
             "of entries with matching values in this field will be merged "
             "across all provided relational databases")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if (args.fields and not args.map_files) or \
        (args.map_files and not args.fields):
        parser.error("error: -m/--mapping and -f/--fields must be supplied "
                     "together")

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

    sample_names = args.names if args.names else [os.path.basename(i) for i \
                   in args.csvs]
    feature_names = args.features
    fields = args.fields

    # Store feature abundances in a dictionary
    totals = 0
    abundances = {}
    for position, csv_file in enumerate(args.csvs):
        with open_io(csv_file, mode='rb') as in_h:
            for row in in_h:
                row = row.decode('utf-8')

                if row.startswith('#'):  #skip comments
                    continue

                row = row.split('\t')

                try:
                    name, abundance = row
                except:
                    raise FormatError("input file {} does not have the "
                                      "correct number of columns. Please "
                                      "verify that that the file is formatted "
                                      "correctly".format(csv_file))

                # Ignore features not found in list of features to include
                if feature_names:
                    if name not in feature_names:
                        continue

                try:
                    abundances[name][position] = float(abundance)
                except KeyError:
                    abundances[name] = [float(0) for i in sample_names]

    if args.map_files:
        # Create header from sample names and desired fields
        header_columns = fields + sample_names

        # Read in the relational databases
        mapping = {}
        all_fields = fields + [args.dup_field] if args.dup_field else fields
        for map_file in args.map_files:
            json_map = json.load(open_io(map_file))

            # Iterate over JSON fiel
            for item in json_map:
                entry = {k: json_map[item][k] for k in json_map[item].keys() & \
                        set(all_fields)}
                mapping[item] = entry

        if args.dup_file:
            mapping = merge_dbs.derep_by_file(mapping, args.dup_file)
        elif args.dup_field:
            mapping = merge_dbs.derep_by_field(mapping, args.dup_field)

    else:
        header_columns = sample_names
        mapping = None

    # Output header
    header = "Feature\t{}\n".format('\t'.join(header_columns))
    out_h(header.encode('utf-8'))

    # Output feature abundances by sample
    list_type = type(list())
    for feature in sorted(abundances):
        entries = []
        if mapping:
            for field in fields:
                try:
                    entry = mapping[feature][field]
                except KeyError:
                    entry = "NA"
                except AttributeError:
                    entry = "NA"
                if not entry:
                    entry = "NA"

                if type(entry) == list_type:
                    entry = ';'.join(entry)

                entries.append(entry)

        entries = entries + ['{0:.2f}'.format(i) for i in abundances[feature]]
        out_h("{}\t{}\n".format(feature, '\t'.join(entries)).encode('utf-8'))

    # Output statistics
    f_totals = len(abundances)
    s_totals = len(sample_names)
    print("Features processed:", file=sys.stderr)
    print("  - feature totals:\t{!s}".format(f_totals), file=sys.stderr)
    print("  - samples merged:\t{!s}".format(s_totals), file=sys.stderr)
    print("", file=sys.stderr)

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to merge {!s} samples"\
          .format(total_time, s_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
