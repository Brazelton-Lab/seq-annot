#! /usr/bin/env python
"""
Annotate a GFF3 file of putative protein-coding genes using the results of a 
homology search to a database of genes with known or predicted function.

Usage:
    annotate_features [options] [-m in.json] [-o out.gff] -b in.b6 in.gff

Required input is a GFF3 file of predicted protein-coding genes and a tabular 
file of pairwise alignments (B6 format). Optional input is a relational 
database in JSON format containing additional information on a hit. It is 
assumed that the query IDs in the alignment file are the sequence IDs in the 
GFF3 file combined with the second part of the ID tag, separated by an 
underscore.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension to 
the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to 
standard output (stdout).

Copyright:

    annotate_features annotates predicted features based on homology
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

from __future__ import print_function

from arandomness.argparse import Open, ParseSeparator
import argparse
from bio_utils.iterators import b6_iter, gff3_iter
import json
from seq_annot.seqio import open_input
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.0.1"


def do_nothing(*args):
    pass


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('gff',
        metavar='in.gff',
        action=Open,
        mode='rb',
        default=sys.stdin,
        help="input feature predictions in GFF3 format")
    parser.add_argument('-b', '--b6',
        action=ParseSeparator,
        sep=',',
        help="input one or more best-hit alignment files in B6/M8 format. The "
             "alignment format should be provided to -s/--specifiers if other "
             "than the default BLAST+ tabular format. If a given feature has "
             "a match in more than one alignment file, precedence will be "
             "determined by input order")
    parser.add_argument('-m', '--mapping',
        metavar='in.json',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input one or more relational databases in JSON format containing "
        "attribute fields")
    parser.add_argument('-o', '--out',
        metavar='out.gff',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output annotated features in GFF3 format [default: output to "
             "stdout]")
    parser.add_argument('-s', '--specifiers',
        dest="format",
        action=ParseSeparator,
        sep=",",
        help="input alignment format, ordered by field specifier [default: "
             "qaccver, saccver, pident, length, mismatch, gapopen, qstart, "
             "qend, sstart, send, evalue, bitscore]")
    parser.add_argument('-f', '--fields',
        metavar='FIELD [,FIELD,...]',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the relational database "
             "that should be included as attributes, if relevant")
    parser.add_argument('-p', '--precedence',
        choices=['order', 'quality'],
        default='quality',
        help="method to resolve conflicting annotations [default: quality]. "
             "Options are 'order' and 'quality'. If 'order', precedence will "
             "be determined by input file order. If 'quality', precedence "
             "will be given to the match with the highest quality alignment "
             "as determined by bitscore")
    parser.add_argument('-d', '--default',
        metavar='STR',
        dest='prod_def',
        help="default product for features without associated product "
             "information")
    output_control = parser.add_argument_group(title="output control options")
    output_control.add_argument('--filter',
        action='store_true',
        help="do not output unannotated features [default: no filtering]")
    output_control.add_argument('--keep-attrs',
        dest='keep',
        action='store_true',
        help="do not overwrite existing attributes, but append new attributes "
             "to the end of the list [default: overwrite existing]")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if not (args.fields and args.map_files):
        parser.error("error: -m/--mapping and -f/--fields must be supplied "
                     "together")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('annotate_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    out_h = args.out.write

    if args.map_files:
        mapping = {}
        for map_file in args.map_files:
            json_map = json.load(map_file)
            mapping = {**json_map, **mapping}
    else:
        mapping = None

    map_fields = args.fields
    match_precedence = args.precedence
    specifiers = args.format
    feature_type = args.ftype
    default_product = args.def_prod

    # Initiate statistics variables
    hits_totals = 0
    conflict_totals = 0
    no_map = 0
    no_fields = {}
    for map_field in map_fields:
        map_fields[map_field] = 0

    # Parse results of the homology search
    hits = {}  #store matches and additional attributes
    for b6 in args.b6:

        with open_input(b6) as b6_h:
            for hit in b6_iter(b6_h, header=specifiers):

                hits_totals += 1

                subject = hit.subject

                add_annots = []
                if mapping:
                    try:
                        sub_entry = mapping[subject]
                    except KeyError:
                        no_map += 1
                    else:
                        for field in map_fields:
                            try:
                                add_annots.append((field, sub_entry[field]))
                            except KeyError:
                                no_fields[field] += 1

                query = hit.query

                if query in hits:
                    conflict_totals += 1

                    if match_precedence == 'order':
                        # Preserve only best hit on first come basis
                        continue
                    else:
                        # Determine which match has the best alignment score
                        prev_score = hits[query]['bitscore']

                        if prev_score >= hit.bitscore:
                            continue  #existing match is best
                        else:
                            hits[query] = {'bitscore': hit.bitscore, 'attr': \
                                           [('Name', name)] + add_annots}

                else:
                    hits[query] = {'bitscore': hit.bitscore, 'attr': [('Name', \
                                   name)] + add_annots}

    # Parse GFF file, adding annotations from homology search and relational db
    gff_totals = 0
    passed_totals = 0
    feature_type_totals = 0
    for entry in gff3_iter(args.gff, parse_attr=True, headers=True):
        try:
            seq_id = entry.seqid
        except AttributeError:
            if entry.startswith('##'):
                out_h("{}\n".format(entry))
                continue
            else:
                continue  #don't output comments

        gff_totals += 1

        feature_id = entry.attributes['ID'].split('_')[-1]  #id is second value

        if not args.keep:
            entry.attributes.clear()

        # Skip features of other type
        if feature_type:
            if entry.type == feature_type:
                feature_type_totals += 1
            else:
                unique_id = "f{!s}".format(passed_totals)
                entry.attributes['ID'] = unique_id

                out_h(entry.write())
                continue

        try:
            attrs = hits["{}_{}".format(seq_id, feature_id)]['attr']
        except KeyError:
            if args.filter:
                continue
        else:
            for attr in attrs:
                entry.attributes[attr[0]] = attr[1]  #(name, value)

            # Add default product info if not already available
            if 'product' not in entry.attributes and default_product:
                entry.attributes['product'] = default_product

        passed_totals += 1

        unique_id = "f{!s}".format(passed_totals)
        entry.attributes['ID'] = unique_id

        out_h(entry.write())

    # Calculate and print statistics
    if no_map > 0:
        print("warning: there were {!s} alignments that did not have a "
              "corresponding mapping entry\n".format(no_map), file=sys.stderr)

    for map_field in no_fields:
        no_map_size = no_fields[map_field]
        if no_map_size > 0:
            print("warning: there were {!s} alignments that did not have an "
                  "entry with field {} in a mapping file\n"\
                  .format(no_map_size, field), file=sys.stderr)

    print("Total alignments processed:\t{!s}".format(hits_totals), \
          file=sys.stderr)
    print("Total features processed:\t{!s}".format(gff_totals), file=sys.stderr)
    if feature_type:
        print("  - features of relevant type:\t{!s}"\
              .format(feature_type_totals), file=sys.stderr)
    print("  - features aligned to a single reference:\t{!s}"\
          .format(passed_totals), file=sys.stderr)
    print("  - features aligned to more than one reference:\t{!s}\n"\
          .format(conflict_totals), file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to annotate {!s} features\n"\
          .format(total_time, gff_totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
