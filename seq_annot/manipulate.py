#! /usr/bin/env python
"""
Manipulate abundance table. Can add database fields or filter by cell value.

Usage:
    manipulate_abund in.csv

Required input is one or more CSV files containing feature abundances.

The compression algorithm is automatically detected for input files based on
the file extension. To compress output, add the appropriate file extension
to the output file name (e.g. .gz, .bz2).

Copyright:

    manipulate_abund manipulate an abundance table
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
from seq_annot.seqio import open_io
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = '0.1.0'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('abunds',
        metavar='in.csv',
        action=Open,
        mode='rb',
        help="input tab-separated feature abundance table")
    parser.add_argument('-m', '--mapping',
        metavar='in.json',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input one or more relational databases in JSON format")
    parser.add_argument('-o', '--out',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output modified feature abundance table [default: output to "
             "stdout]")
    parser.add_argument('-f', '--fields',
        metavar='FIELD [,FIELD,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database "
             "that should be added to or removed from the abundance table")
    parser.add_argument('-r', '--filter',
        metavar='FIELD:VALUE [FIELD:VALUE ...]',
        nargs='*',
        help="list of field-value pairs. Features containing the field values "
             "will be filtered out of the abundance table. Arguments must be "
             "quoted if spaces are contained in either the field or value")
    args = parser.parse_args()

    if (args.map_files and not args.fields) or \
        (args.fields and not args.map_files):
        parser.error("error: -m/--mapping and -f/--fields must be supplied "
                     "together")

    # Speedup tricks
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
    fields = args.fields if args.fields else []

    filts = {j[0]: j[1] for j in [i.split(':', 1) for i in args.filter]} if \
            args.filter else []

    # Load databases
    mapping = {}
    if args.map_files:
        for map_file in args.map_files:
            json_map = json.load(open_io(map_file))
            for item in json_map:
                entry = json_map[item]
                entry = {k: entry[k] for k in entry.keys() & set(fields)}

                mapping[item] = entry

    header = args.abunds.readline().decode()
    header = header.strip().split('\t')

    # Insert new fields into table, if applicable
    header = header[0:1] + fields + header[1:]

    # Find column number of fields to filter
    indexes = {}
    for filt in filts:
        if filt not in header:
            print("warning: unknown column name '{}' provided to --filter"\
                  .format(filt), file=sys.stderr)
        else:
            index = header.index(filt)
            try:
                indexes[index].append(filts[filt])
            except KeyError:
                indexes[index] = [filts[filt]]

    feature_ind = header.index('Feature')

    # Output header
    out_h.write('{}\n'.format('\t'.join(header)))

    # Read abundance table
    for line in args.abunds:
        line = line.decode()
        if line.startswith('#'):
            continue
        
        split_line = line.strip().split('\t')
        feature = split_line[feature_ind]

        #Insert new field values
        new_values = []
        for field in fields:
            try:
                field_val = mapping[feature][field]
            except KeyError:
                field_val = 'NA'

            if type(field_val) == list_type:
                field_val = ';'.join(field_val)

            if not field_val:
                field_val = 'NA'

            new_values.append(field_val)
        
        split_line = split_line[0:1] + new_values + split_line[1:]

        # Filter table
        fail = 0
        for index in indexes:
            value = split_line[index].split(';')
            if set(value).intersection(indexes[index]):
                fail = 1
                break
        if fail:
            continue
        else:
            try:
                out_h.write('{}\n'.format('\t'.join(split_line)))
            except TypeError:
                print(split_line)
                sys.exit(1)

    # Calculate and print program run-time info
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to modify the abundance table"\
          .format(total_time), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
