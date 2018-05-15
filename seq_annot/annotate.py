#! /usr/bin/env python
"""
Annotate a GFF3 file of putative protein-coding genes using the results of a 
homology search to a database of genes with known or predicted function.

Required input is a GFF3 file of predicted protein-coding genes and a tabular 
file of pairwise alignments (B6 format). Optional input is a relational 
database in JSON format containing additional information on a hit. It is 
assumed that the query IDs in the alignment file are the sequence IDs in the 
GFF3 file combined with the second part of the ID tag, separated by an 
underscore.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension to 
the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to 
standard output (stdout). Input is taken from standard input (sdtin) by default.

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
__version__ = "0.2.1"


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
        help="input feature predictions in GFF3 format [default: input from "
             "stdin]")
    parser.add_argument('-b', '--b6',
        metavar='in.b6 [,in.b6,...]',
        action=ParseSeparator,
        sep=',',
        help="input one or more best-hit alignment files in B6/M8 format. The "
             "alignment format should be provided to -s/--specifiers if other "
             "than the default BLAST+ tabular format. If a given feature has "
             "a match in more than one alignment file, precedence will be "
             "determined by input order")
    parser.add_argument('-m', '--mapping',
        metavar='in.json [,in.json,...]',
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
    parser.add_argument('-c', '--conflict',
        choices=['order', 'quality'],
        default='quality',
        help="method to resolve conflicting annotations [default: quality]. "
             "Options are 'order' and 'quality'. If 'order', precedence will "
             "be determined by input file order. If 'quality', precedence "
             "will be given to the match with the highest quality alignment, "
             "as determined by bitscore")
    parser.add_argument('-t', '--type',
        dest='ftype',
        help="feature type (3rd column in GFF file) to annotate [default: "
             "annotate features of all type]. All features of other type will "
             "be ignored")
    parser.add_argument('-p', '--product',
        metavar='STR',
        dest='prod_def',
        help="default product for features without associated product "
             "information in relational database")
    output_control = parser.add_argument_group(title="output control options")
    output_control.add_argument('-d', '--discarded',
        metavar='out.log',
        action=Open,
        mode='wt',
        help="output information on discarded annotations")
    output_control.add_argument('--filter',
        action='store_true',
        help="do not output unannotated features [default: no filtering]")
    output_control.add_argument('--clear',
        dest='clear_attrs',
        action='store_true',
        help="overwrite existing attributes [default: append new attributes "
             "to the end of existing attributes]")
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
    print("{} {!s}".format('annotate_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    out_h = args.out.write
    out_log = args.log.write if args.log else do_nothing
    out_log("#Kept\tDiscarded\tReason\n")

    if args.map_files:
        mapping = {}
        for map_file in args.map_files:
            json_map = json.load(open_input(map_file))
            mapping = {**json_map, **mapping}
    else:
        mapping = None

    map_fields = args.fields
    match_precedence = args.precedence
    specifiers = args.format
    default_product = args.prod_def
    feature_type = args.ftype

    # Initiate statistics variables
    hits_totals = 0
    conflict_totals = 0
    no_map = 0
    no_fields = {}
    for map_field in map_fields:
        no_fields[map_field] = 0

    # Parse results of the homology search
    hits = {}  #store matches and additional attributes
    for b6 in args.b6:

        with open_input(b6) as b6_h:
            for hit in b6_iter(b6_h, header=specifiers):

                hits_totals += 1

                query = hit.query

                if query in hits:
                    conflict_totals += 1

                    if match_precedence == 'order':
                        # Preserve only best hit on first come basis
                        out_log("{}\t{}\t{}\n".format(hits[query].subject, \
                                hit.subject, 'order'))
                        continue
                    else:
                        # Determine which match has the best alignment score
                        prev_score = hits[query].bit_score

                        if prev_score >= hit.bit_score:
                            out_log("{}\t{}\t{}\n".format(hits[query].subject, \
                                    hit.subject, 'score'))
                            continue  #existing match is best
                        else:
                            out_log("{}\t{}\t{}\n".format(hit.subject, \
                                    hits[query].subject, 'score'))
                            hits[query] = hit

                else:
                    hits[query] = hit

    # Parse GFF file, adding annotations from homology search and relational db
    gff_totals = 0
    annotated_totals = 0
    feature_type_totals = 0
    unknown_totals = 0
    prev_name = None
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

        # clear existing attributes if directed
        if args.clear_attrs:
            entry.attributes.clear()

        # Keep track of output ID attribute components
        if seq_id == prev_name:
            seq_count += 1
        else:
            seq_count = 1
            prev_name = seq_id

        # Provide new ID for the feature
        unique_id = "{}_{!s}".format(seq_id, seq_count)
        entry.attributes['ID'] = unique_id

        # Annotate features of a given type only
        if feature_type:
            if feature_type != entry.type:
                if not args.filter:
                    out_h(entry.write())
                continue
            else:
                feature_type_totals += 1

        # Annotate features using attributes field
        try:
            # Query name should be sequenceID_featureID
            hit = hits["{}_{}".format(seq_id, feature_id)]
        except KeyError:
            # Discard feature if unable to annotate
            if args.filter:
                continue
        else:
            subject = hit.subject

            attrs = [('Name', subject)]
            if mapping:
                try:
                    sub_entry = mapping[subject]
                except KeyError:
                    no_map += 1
                else:
                    for field in map_fields:
                        try:
                            entry_value = sub_entry[field]
                        except KeyError:
                            no_fields[field] += 1
                        else:
                            if entry_value:
                                attrs.append((field, entry_value))

            annotated_totals += 1
            for attr in attrs:
                entry.attributes[attr[0]] = attr[1]  #(name, value)

        # Add default to product attribute if doesn't exist
        if 'product' not in entry.attributes and default_product:
            entry.attributes['product'] = default_product

        # Write features to new GFF3 file
        out_h(entry.write())

    # Calculate and print statistics
    if no_map > 0:
        print("warning: there were {!s} alignments that did not have a "
              "corresponding mapping entry\n".format(no_map), file=sys.stderr)

    for map_field in no_fields:
        no_map_size = no_fields[map_field]
        if no_map_size > 0:
            print("warning: there were {!s} alignments that did not have an "
                  "entry with field '{}' in a mapping file\n"\
                  .format(no_map_size, field), file=sys.stderr)

    print("Alignments processed:\t{!s}".format(hits_totals), \
          file=sys.stderr)
    print("Features processed:\t{!s}".format(gff_totals), file=sys.stderr)
    if feature_type:
        print("  - of relevant type:\t{!s}"\
              .format(feature_type_totals), file=sys.stderr)
    print("  - single aligned reference:\t{!s}"\
          .format(annotated_totals), file=sys.stderr)
    print("  - multiple aligned references:\t{!s}\n"\
          .format(conflict_totals), file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to annotate {!s} features\n"\
          .format(total_time, gff_totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
