#! /usr/bin/env python
"""
Annotate a GFF3 file of putative protein-coding genes using the results of a 
homology search to a database of genes with known or predicted function.

Usage:
    annotate_features [options] [-m in.json] [-o out.gff] -b in.b6 in.gff

Required input is a GFF3 file of predicted protein-coding genes and a tabular 
file of pairwise alignments (B6 format). Optional input is a relational 
database in JSON format containing additional information on a hit.

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
        action=Open,
        mode='rb',
        help="input relational database in JSON format containing attribute "
        "fields")
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
        help="fields from the relational database that should be included as "
             "attributes, if relevant")
    parser.add_argument('-p', '--precedence',
        choices=['order', 'quality'],
        default='quality',
        help="method to resolve conflicting annotations [default: quality]. "
             "Options are 'order' and 'quality'. If 'order', precedence will "
             "be determined by input file order. If 'quality', precedence "
             "will be given to the match with the highest quality alignment "
             "as determined by e-value")
    parser.add_argument('--filter',
        action='store_true',
        help="do not output unannotated features [default: no filtering]")
    parser.add_argument('--keep-attrs',
        dest='keep',
        action='store_true',
        help="do not overwrite existing attributes, but append new attributes "
             "to the end of the list [default: overwrite existing]")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()


    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('annotate_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)


    # Track program run-time
    start_time = time()


    # Assign variables based on user inputs
    map_fields = args.fields
    match_precedence = args.precedence
    specifiers = args.format

    if args.mapping:
        mapping = json.load(args.mapping)


    hits_totals = 0
    conflict_totals = 0
    # Parse results of the homology search
    hits = {}  #store matches and additional attributes
    for b6 in args.b6:

        with open_input(b6) as b6_h:
            for hit in b6_iter(b6_h, header=specifiers):

                hits_totals += 1

                name = hit.subject
                add_annots = []

                if mapping:
                    try:
                        sub_entry = mapping[name]
                    except KeyError:
                        print("error: {} not found in the relational "
                              "database".format(name), file=sys.stderr)
                        sys.exit(1)

                    for i in map_fields:
                        try:
                            add_annots.append((i, sub_entry[i]))
                        except KeyError:
                            print("error: field {} not found in the "
                                  "relational database", file=sys.stderr)
                            sys.exit(1)

                if hit.query in hits:
                    conflict_totals += 1

                    if match_precedence == 'order':
                        # Preserve only best hit on first come basis
                        continue
                    else:
                        # Determine which match has the best alignment score
                        prev_score = hits[hit.query]['evalue']

                        if prev_score <= hit.evalue:
                            continue  #existing match is best
                        else:
                            hits[hit.query] = {'evalue': hit.evalue, 'attr':\
                                               [('Name', name)] + add_annots}

                else:
                    hits[hit.query] = {'evalue': hit.evalue, 'attr': [('Name',\
                                       name)] + add_annots}


    # Parse GFF file, adding annotations from homology search and relational db
    gff_totals = 0
    passed_totals = 0
    for entry in gff3_iter(args.gff, parse_attr=True, headers=True):
        try:
            seq_id = entry.seqid
        except AttributeError:
            if entry.startswith('##'):
                args.out.write("{}\n".format(entry))
                continue
            else:
                continue  #don't output comments

        gff_totals += 1

        feature_id = entry.attributes['ID'].split('_')[-1]  #id is second value

        if not args.keep:
            entry.attributes.clear()
            entry.attributes['ID'] = ''

        try:
            attrs = hits["{}_{}".format(seq_id, feature_id)]['attr']
        except KeyError:
            if args.filter:
                continue
            else:
                entry.attributes['product'] = "hypothetical protein"
        else:
            for attr in attrs:
                entry.attributes[attr[0]] = attr[1]

        if 'product' not in entry.attributes:
            entry.attributes['product'] = "hypothetical protein"

        passed_totals += 1

        unique_id = "f{!s}".format(passed_totals)
        entry.attributes['ID'] = unique_id

        args.out.write(entry.write())


    # Calculate and print statistics
    print("Total matches processed:\t{!s}".format(hits_totals), file=sys.stderr)
    print("Total features processed:\t{!s}".format(gff_totals), file=sys.stderr)
    print("  - features matching a single reference:\t{!s}".format(passed_totals), \
          file=sys.stderr)
    print("  - features with more than one matching reference:\t{!s}\n"\
          .format(conflict_totals), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to annotate {!s} features\n"\
          .format(total_time, gff_totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
