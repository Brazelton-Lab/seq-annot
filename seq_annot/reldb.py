#! /usr/bin/env python
"""
Perform various manipulations on one or more relational databases, including: 
entry merging, database subsetting, and format conversion.

Usage:
    reldbs in.json [in.json ...]

Required input is one or more JSON-formatted database files.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2).

Copyright:

    reldbs manipulate relational databases and database entries
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
import json
import os
from seq_annot.db import *
from seq_annot.seqio import open_io, write_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.4.0'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('map_files',
        metavar='in.json [in.json ...]',
        nargs='+',
        help="input one or more relational databases in JSON format")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output modified relational database [default: output to stdout]")
    merge_group = parser.add_argument_group("merging arguments")
    derep_group = merge_group.add_mutually_exclusive_group()
    derep_group.add_argument('-c', '--dup-file',
        metavar='in.tsv',
        dest='dup_file',
        action=Open,
        mode='rt',
        help="input tab-separated file of database entries to merge. The "
             "input relational databases should be provided as columns "
             "and entries to merge between databases as rows. Lines starting "
             "with '#' will be ignored.")
    derep_group.add_argument('-d', '--dup-field',
        metavar='FIELD',
        dest='dup_field',
        help="database field by which entries should be merged. The content "
             "of the entries with matching values in this field will be "
             "combined across all relational databases input")
    subset_group = parser.add_argument_group("subsetting arguments")
    subset_group.add_argument('-i', '--ids',
        metavar='ID [,ID,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of database entry IDs. Only IDs "
             "in the list will be output [default: output all]")
    subset_group.add_argument('-s', '--subset',
        metavar='FIELD:PATTERN',
        dest='subset_terms',
        action='append',
        help="pattern matching criteria used to extract entries from the "
             "database. Can provide multiple criteria through repeated use "
             "of the argument. The field tag is distinguished from the search "
             "pattern by the first colon character encountered. Entries with "
             "field value matching the pattern will be included in the output")
    subset_group.add_argument('-r', '--filter',
        metavar='FIELD:PATTERN',
        dest='filter_terms',
        action='append',
        help="pattern matching criteria used to filter entries from the "
             "database. Can provide multiple criteria through repeated use "
             "of the argument. The field tag is distinguished from the search "
             "pattern by the first colon character encountered. Features with "
             "field value matching the pattern will be excluded from the "
             "output")
    subset_group.add_argument('--case',
        dest='match_cs',
        action='store_true',
        help="pattern matching should be case-sensitive [default: "
             "pattern matching is case-insensitive]")
    output_group = parser.add_argument_group("output control arguments")
    output_group.add_argument('-f', '--fields',
        metavar='FIELD [,FIELD,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database to "
             "be included in the output")
    output_group.add_argument('--csv',
        action='store_true',
        help="output database as tabular CSV")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('reldbs', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user input
    out_h = args.out
    as_csv = args.csv

    ids = args.ids

    cs = True if args.match_cs else False

    fields = args.fields if args.fields else []
    merge_field = [args.dup_field] if args.dup_field else []
    subset_fields = [i.split(':', 1)[0] for i in args.subset_terms] \
                    if args.subset_terms else []
    filter_fields = [i.split(':', 1)[0] for i in args.filter_terms] \
                    if args.filter_terms else []

    all_fields = fields + subset_fields + filter_fields + merge_field

    merge = True if (args.dup_field or args.dup_file) else False

    if args.dup_file:
        derep_algo = derep_by_file
        by = args.dup_file
    else:
        derep_algo = derep_by_field
        by = args.dup_field

    # Load databases
    mapping = load_dbs(args.map_files, fields=all_fields)
    nentry = len(mapping)

    # Keep only entries in list of IDs supplied by user
    if ids:
        unwanted = set(mapping.keys()) - set(ids)
        for unwanted_key in unwanted:
            del mapping[unwanted_key]

        if not len(mapping) > 0:
            print("error: database contains no entries matching the list of "
                "of supplied IDs", file=sys.stderr)
            sys.exit(1)

    # Remove entries from database if values do not match any of the subset
    # patterns
    if args.subset_terms:
        mapping = filter_dbs(mapping, patterns=args.subset_terms, 
            subset=True, case=cs)

    # Remove entries from database if values match any of the filter patterns
    if args.filter_terms:
        mapping = filter_dbs(mapping, patterns=args.filter_terms,
                             subset=False, case=cs)

    # Merge database entries by field or file
    if merge:
        mapping = derep_algo(mapping, by)

    # Output modified database
    if as_csv:
        # Output CSV header
        output_fields = sorted(set(all_fields) - set(merge_field)) if merge_field else all_fields
        header = "#ID\t{}\n".format("\t".join(output_fields))
        write_io(out_h, header)

        for entry in mapping:
            write_io(out_h, "{}\n".format(entry_as_csv(entry, mapping[entry])))
    else:
        json.dump(mapping, out_h, sort_keys=True, indent=4, \
            separators=(',', ': '))

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to parse {!s} entries"\
          .format(total_time, nentry), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
