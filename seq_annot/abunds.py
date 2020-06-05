#! /usr/bin/env python
"""
Combine abundance values of genomic features based on categories defined in an 
annotated GFF or CSV file.

Required input is a tabular file of feature abundances and either a GFF3 file 
of annotated protein-coding genes or a tabular file relating features to some 
higher order organization.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension to 
the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to 
standard output (stdout). Standard input (stdin) can be redirected to one of 
the positional argument if supplied with '-'.

Copyright:

    combine_abunds combines abundances based on defined categories
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

logging.basicConfig(level=logging.ERROR)

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.2.4"

def group_from_gff(infile, id_tag='ID', attr_tag='Alias'):
    """Obtain grouping information from GFF file
    """
    groupings = {}
    gff_reader = GFF3Reader(infile)
    for entry in gff_reader.iterate(parse_attr=True):
        try:
            feature_id = entry.attributes[id_tag]
        except KeyError:
            raise FormatError("GFF entry for region {} does not have ID tag "\
                "{} required to link back to the abundance file"\
                .format(entry.seqid, id_tag), file=sys.stderr)

        try:
            category = entry.attributes[attr_tag]
        except KeyError:
            logging.warning("GFF entry for region {} does not have attribute "
                "tag '{}'".format(entry.seqid, attr_tag))
            category = "UNKNOWN"

        groupings[feature_id] = category

    return(groupings)

def group_from_csv(infile, *args):
    """Obtain grouping information from CSV file
    """
    groupings = {}
    for nline, row in enumerate(infile):
        row = row.decode('utf-8')

        if row.startswith('#'):  #skip comments
            continue

        row = row.rstrip().split('\t')

        try:
            feature_id, category = row[0], row[1]
        except ValueError:
            raise FormatError("{}, line {!s}: Incorrect number of "
                "columns provided. Please verify file format"\
                .format(infile.name, nline))

        groupings[feature_id] = category

    return(groupings)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('in_abund',
        metavar='in.abund',
        help="input tabular file of feature abundances. Use '-' to indicate "
            "that input should be taken from standard input (stdin)")
    ex_group = parser.add_mutually_exclusive_group()
    ex_group.add_argument('-g', '--gff',
        metavar='in.gff',
        dest='in_gff',
        action=Open,
        mode='rb',
        help="input feature annotations in GFF3 format")
    ex_group.add_argument('-t', '--tsv',
        metavar='in.csv',
        dest='in_csv',
        action=Open,
        mode='rb',
        help="input tabular file grouping feature into desired categories. The "
            "first column should contain feature IDs that correspond to the "
            "IDs in the abundance file. The second column should contain "
            "feature groupings")
    parser.add_argument('-o', '--out',
        metavar='out.abund',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output tabular file of combined abundances [default: output "
            "to stdout]")
    group1 = parser.add_argument_group("GFF options")
    group1.add_argument('-a', '--attr',
        metavar='TAG',
        dest='attr',
        help="attribute tag used to combine feature abundances [default: "
            "None]. If left off, abundances will be combined for a given "
            "genomic region")
    group1.add_argument('-i', '--id',
        metavar='TAG',
        dest='id',
        default='ID',
        help="attribute tag used to match GFF3 entries to feature IDs in the "
            "abundance file [default: ID]")
    group2 = parser.add_argument_group("Combining options")
    group2.add_argument('-s', '--substring',
        metavar='LEN',
        dest='sub',
        type=int,
        help="matching term is a substring of the feature of length LEN")
    group2_mut = group2.add_mutually_exclusive_group()
    group2_mut.add_argument('-m', '--method',
        metavar='[sum|mean|median|occur]',
        dest='method',
        choices=["sum", "mean", "median", "occur"],
        help="method used to combine feature abundance values [default: sum]")
    group2_mut.add_argument('-c', '--cat',
        metavar='SEP',
        dest='cat',
        type=str,
        help="instead of combining abundances, concatenate feature with "
            "feature category using separator character SEP")
    parser.add_argument('--verbose',
        action='store_true',
        help="increase output verbosity")
    args = parser.parse_args()

    logging.basicConfig(level=logging.ERROR)
    if args.verbose:
        logging.basicConfig(level=logging.WARNING)

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
    if args.in_abund == '-':
        in_abund = sys.stdin
    else:
        in_abund = open_io(args.in_abund, mode='rb')

    out_h = args.out
    attr_tag = args.attr
    id_tag = args.id
    separator = args.cat
    substring = args.sub

    if args.method:
        if args.method == "mean":
            method = mean
        elif args.method == "median":
            method = median
        elif args.method == "occur":
            method = occurrence
        else:
            method = sum

    if args.in_gff:
        groups = group_from_gff(args.in_gff, id_tag, attr_tag)
    else:
        groups = group_from_csv(args.in_csv)

    # Iterate over abundance file
    totals = 0
    ncats = 0
    abunds = {}
    for nline, row in enumerate(in_abund):
        if (isinstance(row, bytes)):
            row = row.decode('utf-8')

        if row.startswith('#'):  #skip comments
            continue
        else:
            totals += 1

        row = row.rstrip().split('\t')

        try:
            name = row[0]
        except ValueError:
            infile = os.path.basename(in_abund.name)
            raise FormatError("{}, line {!s}: Incorrect number of "
                "columns provided. Please verify file format"\
                .format(infile, nline))

        if substring:
            name_lookup = name[0: substring]
        else:
            name_lookup = name

        try:
            category = groups[name_lookup]
        except KeyError:
            logging.warning("no grouping information found for {}".format(name))
            category = "NA"
        else:
            ncats += 1

        if args.cat:
            try:
                cols = "\t".join(row[1:])
            except ValueError:
                infile = os.path.basename(in_abund.name)
                raise FormatError("{}, line {!s}: Incorrect number of "
                    "columns provided. Please verify file format"\
                    .format(infile, nline))

            output = "{}{}{}\t{}\n".format(name, separator, category, cols)
            write_io(out_h, output)

        else:
            try:
                value = row[1]
            except ValueError:
                infile = os.path.basename(in_abund.name)
                raise FormatError("{}, line {!s}: Incorrect number of "
                    "columns provided. Please verify file format"\
                    .format(infile, nline))

            try:
                abunds[category].append(float(value))
            except KeyError:
                abunds[category] = [float(value)]
            except ValueError:
                infile = os.path.basename(args.in_abund)
                raise FormatError("{}, line {!s}: abundance value must be "
                    "either an integer or float".format(infile, nline))

    # Output merged abundances
    if not args.cat:
        for category in abunds:
            final_value = method(abunds[category])
            write_io(out_h, "{}\t{!s}\n".format(category, final_value))

    # Output statistics
    nfeatures = len(groups.keys())
    print("", file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - features combined:\t{!s}".format(totals), file=sys.stderr)
    print("  - output categories:\t{!s}".format(ncats), file=sys.stderr)
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
