#! /usr/bin/env python
"""
Manipulate a feature abundance table by inserting new database fields and 
subsetting via screening of the relational database.

Required input is a feature abundance table in TSV format and a relational
database in JSON format.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to
standard output (stdout). The input abundance table will be taken from 
standard input (sdtin) by default.

Copyright:

    manipulate_abund manipulate a feature abundance table
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

from arandomness.argparse import Open, ParseSeparator
import argparse
import os
from seq_annot.reldb import filter_dbs, load_dbs
from seq_annot.seqio import open_io, write_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.2.4'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('abunds',
        metavar='in.csv',
        action=Open,
        mode='rb',
        default=sys.stdin,
        help="input tab-separated feature abundance table")
    parser.add_argument('-m', '--mapping',
        required=True,
        metavar='in.json [,in.json,...]',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input one or more relational databases in JSON format")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output modified feature abundance table [default: output to "
             "stdout]")
    parser.add_argument('-i', '--insert',
        metavar='FIELD [,FIELD,...]',
        dest='input_fields',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database "
             "that should be added to the abundance table")
    parser.add_argument('-s', '--subset',
        metavar='FIELD:PATTERN',
        dest='search_terms',
        action='append',
        help="pattern-matching criteria used to subset the relational "
             "database. Can provide multiple criteria through repeated use "
             "of the argument. Will separate field from the search pattern on "
             "first encounter of the colon character. Features with field "
             "value matching the pattern will be returned.")
    parser.add_argument('--filter',
        dest='filter',
        action='store_true',
        help="discard entries ")
    args = parser.parse_args()

    if not (args.input_fields or args.search_terms):
        parser.error("-m/--mapping must be used with either "
                     "-i/--insert or -f/--filter")

    # Speedup trick
    list_type = type(list())

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('manipulate_abunds', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user input
    out_h = args.out

    add_fields = args.input_fields if args.input_fields else []
    search_fields = [i.split(':', 1)[0] for i in args.search_terms] \
                    if args.search_terms else []

    keep_match = True if not args.filter else False

    # Load databases
    mapping = load_dbs(args.map_files, fields=add_fields + search_fields)

    # Remove entries from database if values do not match any of the patterns
    mapping = filter_dbs(mapping, patterns=args.search_terms, subset=keep_match)

    # Insert new fields into table header, if applicable
    try:
        header = args.abunds.readline().decode('utf-8')
    except UnicodeDecodeError:
        header = args.abunds.readline()

    header = header.strip().split('\t')
    header = header[0:1] + add_fields + header[1:]

    write_io(out_h, '{}\n'.format('\t'.join(header)))

    # Iterate over feature abundance table
    for row in args.abunds:
        try:
            row = row.decode('utf-8')
        except AttributeError:
            pass

        if row.startswith('#'):  # Skip comments
            continue
        
        split_row = row.strip().split('\t')

        feature_name = split_row[0]

        # Retain only those feature abundances with corresponding database 
        # entries after database filtering
        if search_fields:
            if feature_name not in mapping:
                continue

        # Insert new field values
        new_values = []
        for add_field in add_fields:
            try:
                field_val = mapping[feature_name][add_field]
            except KeyError:  # Feature or entry field not in database
                field_val = 'NA'

            if type(field_val) == list_type:  # Attribute has multiple values
                field_val = ';'.join(field_val)

            if not field_val:
                field_val = 'NA'

            new_values.append(field_val)
 
        # Output row
        output = "\t".join(split_row[0:1] + new_values + split_row[1:])
        write_io(out_h, "{}\n".format(output))

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to modify the abundance table"\
          .format(total_time), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
