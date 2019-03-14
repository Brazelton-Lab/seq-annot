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
import json
import os
from seq_annot.reldb import load_dbs
from seq_annot.seqio import open_io, write_io, FormatError
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.4.1'


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
    parser.add_argument('-a', '--alias',
        metavar='FIELD',
        dest='alias_field',
        help="database field containing the alias by which a feature should "
             "be known. Feature abundances will be combined for features with "
             "identical alias")
    parser.add_argument('--header',
        action='store_true',
        help="first line of the input file is a header")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if args.fields and not args.map_files:
        parser.error("error: -m/--mapping and -f/--fields must be supplied "
                     "together")

    if args.alias_field and not args.map_files:
        parser.error("error: -m/--mapping and -a/--alias must be supplied "
                     "together")

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

    sample_names = args.names if args.names else [os.path.basename(i) for i \
                   in args.csvs]
    feature_names = args.features
    fields = args.fields if args.fields else []
    alias_field = args.alias_field
    end_pos = len(fields)

    if args.map_files and (alias_field or fields):
        # Create header from sample names and desired fields
        header_columns = fields + sample_names

        # Read in the relational databases
        all_fields = fields + [args.alias_field] if args.alias_field else fields
        mapping = load_dbs(args.map_files, fields=all_fields)
    else:
        header_columns = sample_names
        mapping = None

    # Store feature abundances in a dictionary
    no_alias = 0
    totals = 0
    abundances = {}
    for position, csv_file in enumerate(args.csvs):
        with open_io(csv_file, mode='rb') as in_h:
            if args.header:
                header = in_h.readline()

            for row in in_h:

                row = row.decode('utf-8')

                if row.startswith('#'):  #skip comments
                    continue
                else:
                    totals += 1

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

                # Search for alias
                try:
                    feature_id = mapping[name][alias_field]
                except KeyError:
                    no_alias += 1
                    feature_id = name
                except (AttributeError, TypeError):
                    feature_id = name

                try:
                    abundances[feature_id][position + end_pos] += float(abundance)
                except KeyError:
                    entries = []
                    for field in fields:
                        try:
                            entry = mapping[name][field]
                        except KeyError:
                            entry = "NA"
                        except AttributeError:
                            entry = "NA"
                        if not entry:
                            entry = "NA"

                        if type(entry) == list_type:
                            entry = ';'.join(entry)

                        entries.append(entry)

                    abundances[feature_id] = entries + [float(0) for i in \
                                                        sample_names]
                    abundances[feature_id][position + end_pos] += float(abundance)

    # Output header
    header = "Feature\t{}\n".format('\t'.join(header_columns))
    write_io(out_h, header)

    # Output feature abundances by sample
    for feature in abundances:
        entries = abundances[feature]
        entries = '\t'.join(entries[0:end_pos] + ['{:g}'.format(i) for i in \
                                                  entries[end_pos:]])
        write_io(out_h, "{}\t{!s}\n".format(feature, entries))

    # Output statistics
    if alias_field:
        print("", file=sys.stderr)
        print("warning: unable to find an alias for {!s} features"\
              .format(no_alias), file=sys.stderr)

    f_totals = len(abundances)
    s_totals = len(sample_names)
    print("", file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - samples merged:\t{!s}".format(s_totals), file=sys.stderr)
    print("  - feature totals:\t{!s}".format(f_totals), file=sys.stderr)
    if alias_field:
        print("  - features combined:\t{!s}".format(totals - f_totals), \
              file=sys.stderr)

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
