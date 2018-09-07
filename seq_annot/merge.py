#! /usr/bin/env python
"""
Merge relplicate entries among relational databases.

Usage:
    merge_reldbs in.json [in.json ...]

Required input is one or more JSON-formatted files of feature abundances.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2).

Copyright:

    merge_reldbs combine database entries
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
from seq_annot.reldb import load_dbs, derep_by_field, derep_by_file
from seq_annot.seqio import open_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.3.0'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('map_files',
        metavar='in.json [in.json ...]',
        nargs='+',
        help="input one or more relational databases in JSON format to merge")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output merged relational database [default: output to stdout]")
    parser.add_argument('-f', '--fields',
        metavar='FIELD [,FIELD,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database "
            "that should be added to the table")
    derep_group = parser.add_mutually_exclusive_group(required=True)
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
        help="database field through which entries should be combined. The "
             "content of entries with matching values in this field will be "
             "merged across all provided relational databases")
    args = parser.parse_args()

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('merge_reldbs', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user input
    out_h = args.out

    fields = args.fields
    if fields and args.dup_field:
        fields.append(args.dup_field)

    # Load databases
    mapping = load_dbs(args.map_files, fields=fields)

    if args.dup_file:
        derep_algo = derep_by_file
        by = args.dup_file
    else:
        derep_algo = derep_by_field
        by = args.dup_field

    mapping = derep_algo(mapping, by)

    json.dump(mapping, out_h, sort_keys=True, indent=4, \
        separators=(',', ': '))

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to merge {!s} entries"\
          .format(total_time, len(mapping)), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
