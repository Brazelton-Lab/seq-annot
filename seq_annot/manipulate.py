#! /usr/bin/env python
"""
Manipulate feature abundance table. Objects can be filtered by searching a 
relational database for specified patterns. The database can also be used to 
merge objects based on shared field values. If desired, an accompanying 
table of contextual information in TSV format will be created from the 
database fields.

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
from seq_annot.reldb import *
from seq_annot.seqio import open_io, write_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.4.2'


def combine_rows(r1: list, r2: list):
    """Perform vector addition for two rows in an abundance table.
    """
    num_type = type(float())

    if len(r1) != len(r2):
        print("error: rows have differing numbers of columns", file=sys.stderr)
        sys.exit(1)

    try:
        merged = list(map(sum, zip(r1, r2)))
    except TypeError:
        print("error: all cells must have numerical values", file=sys.stderr)
        sys.exit(1)

    return merged


def format_entry(mapping: dict, entry_id: str, fields=None):
    """Format relevant database entry for writing
    """
    try:
        entry = mapping[entry_id]
    except KeyError:
        values = ["NA"] * len(fields)
        entry_csv = "{}\t{}\n".format(entry_id, '\t'.join(values))
    else:
        entry_csv = "{}\t{}\n".format(entry_id, entry_as_csv(entry, fields))

    return entry_csv


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
        metavar='out.csv',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output modified feature abundance table [default: output to "
             "stdout]")
    parser.add_argument('-d', '--db',
        metavar='out.csv',
        action=Open,
        mode='wb',
        help="output accompanying contextual information from relational "
             "database. Can be used with -a/--fields to specify which fields"
             "in the database to include in the output")
    parser.add_argument('-c', '--combine',
        metavar='FIELD',
        dest='merge_field',
        help="combine features with shared values in the database field")
    parser.add_argument('-a', '--fields',
        metavar='FIELD [,FIELD,...]',
        dest='input_fields',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database "
             "that should be added to the accompanying table of contextual "
             "information")
    parser.add_argument('-s', '--subset',
        metavar='FIELD:PATTERN',
        dest='subset_terms',
        action='append',
        help="pattern matching criteria used to subset the relational "
             "database. Can provide multiple criteria through repeated use "
             "of the argument. Will separate field from the search pattern on "
             "first encounter of the colon character. Features with field "
             "value matching the pattern will be returned")
    parser.add_argument('-f', '--filter',
        metavar='FIELD:PATTERN',
        dest='filter_terms',
        action='append',
        help="pattern matching criteria used to filter the relational "
             "database. Can provide multiple criteria through repeated use "
             "of the argument. Will separate field from the search pattern on "
             "first encounter of the colon character. Features with field "
             "value matching the pattern will be excluded")
    parser.add_argument('--case',
        dest='match_cs',
        action='store_true',
        help="pattern matching should be case-sensitive")
    args = parser.parse_args()

    if args.map_files and not (args.input_fields or args.subset_terms or \
                               args.merge_field or args.filter_terms):
        parser.error("-m/--mapping must be used with one or more of "
                     "-c/--combine, -d/--db, -s/--subset, or -f/--filter")

    if args.db and not args.input_fields:
        parser.error("-d/--db should be used -a/--fields")

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
    out_db = args.db

    cs = True if args.match_cs else False

    merge_field = args.merge_field
    add_fields = args.input_fields if args.input_fields else []
    subset_fields = [i.split(':', 1)[0] for i in args.subset_terms] \
                    if args.subset_terms else []
    filter_fields = [i.split(':', 1)[0] for i in args.filter_terms] \
                    if args.filter_terms else []

    all_fields = add_fields + subset_fields + filter_fields + [merge_field]

    # Load databases
    mapping = load_dbs(args.map_files, fields=all_fields) if args.map_files \
              else {}

    # Remove entries from database if values do not match any of the subset 
    # patterns
    if args.subset_terms:
        mapping = filter_dbs(mapping, patterns=args.subset_terms, subset=True, 
                             case=cs)

    # Remove entries from database if values match any of the filter patterns
    if args.filter_terms:
        mapping = filter_dbs(mapping, patterns=args.filter_terms, 
                             subset=False, case=cs)

    # Output headers
    try:
        header = args.abunds.readline().decode('utf-8')
    except AttributeError:
        header = args.abunds.readline()

    write_io(out_h, header)

    if out_db and not merge_field:
        db_header = "{}\t{}\n".format(header.split('\t')[0], 
                                      "\t".join(add_fields))
        write_io(out_db, db_header)

    # Iterate over feature abundance table
    merged_data = {}
    for line_pos, row in enumerate(args.abunds):
        try:
            row = row.decode('utf-8')
        except AttributeError:
            pass

        if row.startswith('#'):  #skip comments
            continue
        
        split_row = row.strip().split('\t')

        feature_name = split_row[0]

        # Retain only those feature abundances with corresponding database 
        # entries after database filtering
        if subset_fields or filter_fields:
            if feature_name not in mapping:
                continue

        # Output row or merge based on field value
        if merge_field != None:
            try:
                field_val = mapping[feature_name][merge_field]
            except KeyError:
                field_val = "NA"

            try:
                counts = list(map(float, split_row[1:]))
            except ValueError:
                print("error: found non-numeric abundance values at line {!s}"\
                      .format(line_pos), file=sys.stderr)
                sys.exit(1)

            try:
                existing = merged_data[field_val]
            except KeyError:
                merged_data[field_val] = counts
            else:
                merged_data[field_val] = combine_rows(counts, existing)
        else:
            output = "\t".join(split_row)
            write_io(out_h, "{}\n".format(output))

            # Write corresponding database entries as CSV
            if out_db:
                write_io(out_db, format_entry(mapping, feature_name, add_fields))

    if merge_field != None:
        if out_db:
            mapping = derep_by_field(mapping, field=merge_field)
            add_fields = [i for i in add_fields if i != merge_field]
            db_header = "{}\t{}\n".format(merge_field, "\t".join(add_fields))
            write_io(out_db, db_header)

        num_type = type(float())

        for new_id in merged_data:
            output = ["{:.9g}".format(i) for i in merged_data[new_id]]
            write_io(out_h, "{!s}\t{}\n".format(new_id, "\t".join(output)))

            # Write corresponding database entries as CSV
            if out_db:
                write_io(out_db, format_entry(mapping, new_id, add_fields))

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to modify the abundance table"\
          .format(total_time), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
