#! /usr/bin/env python
"""
Annotate a GFF3 file of putative protein-coding genes using the results of a 
homology search to a database of genes with known or predicted function.

Required input is a GFF3 file of predicted protein-coding genes and one or 
more pairwise alignments in B6 or HMMER format. Optional input is a relational 
database in JSON format containing additional information on a database entry. 
The dataset IDs in the alignment file must match their corresponding value in 
the ID attribute of the GFF.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension to 
the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output to 
standard output (stdout). Standard input (stdin) can be redirected to the 
positional argument if supplied with '-'.

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

import argparse
from bio_utils.iterators import B6Reader, GFF3Reader
import json
from seq_annot.argparse import *
from seq_annot.db import load_dbs
from seq_annot.seqio import *
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.5.8"


def do_nothing(*args):
    pass

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('in_gff',
        metavar='in.gff',
        help="input feature predictions in GFF3 format. Use '-' to indicate "
             "that input should be taken from standard input (stdin)")
    parser.add_argument('-b', '--hits',
        metavar='in.aln [,in.aln,...]',
        action=ParseSeparator,
        sep=',',
        help="input one or more best-hit alignment files in B6/M8 or HMMER "
            "tblout format. The alignment format should be provided to "
            "-s/--specifiers if other than the default BLAST+ tabular or "
            "HMMER tblout format")
    parser.add_argument('-r', '--format',
        dest="aln_format",
        metavar='FORMAT',
        choices=["b6", "hmm"],
        default="b6",
        help="format of the files containing the homology search results "
            "[default: b6]. Options are b6 or hmm")
    parser.add_argument('-m', '--mapping',
        metavar='in.json [,in.json,...]',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input one or more relational databases in JSON format containing "
            "attribute fields. Multiple input files can be provided by "
            "separating them with a comma and no spaces")
    parser.add_argument('-o', '--out',
        metavar='out.gff',
        action=Open,
        mode='wb',
        default=sys.stdout,
        help="output annotated features in GFF3 format [default: output to "
            "stdout]")
    parser.add_argument('-i', '--id',
        dest='id',
        default='ID',
        help="attribute name used to match GFF3 entries to sequences from the "
            "homology search [default: ID]")
    parser.add_argument('-s', '--specifiers',
        dest="format",
        action=ParseSeparator,
        sep=",",
        default=["qaccver", "saccver", "pident", "length", "mismatch", 
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", 
            "bitscore"], 
        help="input alignment format, ordered by field specifier [default: "
            "qaccver, saccver, pident, length, mismatch, gapopen, qstart, "
            "qend, sstart, send, evalue, bitscore]")
    parser.add_argument('-f', '--fields',
        metavar='FIELD [,FIELD,...]',
        dest='fields',
        action=ParseSeparator,
        sep=',',
        help="comma-separated list of fields from the provided relational "
            "database that should be included as attributes. The field "
            "'Alias' should not be used, as it is reserved for annotation")
    parser.add_argument('-c', '--conflict',
        dest='precedence',
        choices=['order', 'quality'],
        default='quality',
        help="method to resolve conflicting annotations [default: quality]. "
            "Options are 'order' and 'quality'. If 'order', precedence will "
            "be determined by input file order. If 'quality', precedence "
            "will be given to the match with the highest quality alignment, "
            "determined by bitscore")
    parser.add_argument('-t', '--type',
        dest='ftype',
        help="feature type (3rd column in GFF file) to annotate [default: "
            "annotate features of all type]. All features of other type will "
            "be ignored")
    hmm_parser = parser.add_argument_group(title="HMM input options")
    hmm_parser.add_argument('--reverse',
        action='store_true',
        help="dataset sequence assumed the role of target in the homology "
            "search (i.e. hmmsearch was used)")
    hmm_parser.add_argument('--domain',
        action='store_true',
        help="use scores from the single best-scoring domain")
    output_control = parser.add_argument_group(title="output control options")
    output_control.add_argument('-p', '--product',
        metavar='STR',
        dest='prod_def',
        help="default product for features without associated product "
            "information in relational database (e.g. hypothetical protein)")
    output_control.add_argument('-a', '--alias',
        metavar='Field:Alias',
        dest='field_alias',
        action='append',
        help="alias to use in place of the field name when outputting "
            "attributes (e.g. gene:Name). Can provide multiple aliases "
            "through repeated use of the argument")
    output_control.add_argument('-d', '--discarded',
        metavar='out.log',
        dest='log',
        action=Open,
        mode='wb',
        help="output information on discarded annotations")
    output_control.add_argument('--filter',
        action='store_true',
        help="do not output unannotated features [default: no filtering]")
    output_control.add_argument('--clear',
        dest='clear_attrs',
        action='store_true',
        help="overwrite existing attributes [default: append new attributes "
            "to the end of existing attributes]. Attributes with predefined "
            "meaning according to the GFF3 specifications or reserved "
            "attributes (those that begin with an uppercase character) are "
            "exempt")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if args.fields and not args.map_files:
        parser.error("-m/--mapping must be supplied when -f/--fields is used")

    if args.field_alias and not (args.map_files and args.fields):
        parser.error("-m/--mapping and -f/--fields must be supplied when "
            "-a/--alias is used")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('annotate_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Speedup tricks
    string_type = type(str())

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    out_h = args.out
    out_log = args.log.write if args.log else do_nothing
    out_log("#Kept\tDiscarded\tReason\n".encode('utf-8'))

    if args.in_gff == '-':
        in_gff = sys.stdin
    else:
        in_gff = open_io(args.in_gff, mode='rb')

    id_attr = args.id
    aln_format = args.aln_format
    map_fields = [i for i in args.fields if i != "Alias"] \
        if args.fields else []  #Alias is reserved for annotation
    attr_alias = args.field_alias if args.field_alias else []
    match_precedence = args.precedence
    specifiers = args.format
    default_product = args.prod_def
    feature_type = args.ftype

    dom = args.domain
    revq = args.reverse

    # Load mappings, if provided
    mapping = load_dbs(args.map_files) if args.map_files else None

    no_fields = {}  # Mapping-related statistics
    for map_field in map_fields:
        no_fields[map_field] = 0

    # Dictionary of fields and their respective aliases
    attr_names = {}
    alias_fields = [i.split(':')[0] for i in attr_alias]
    all_fields = attr_alias + [i for i in map_fields if i not in alias_fields]
    if default_product and 'product' not in all_fields:
        all_fields.append('product')
    for attr in all_fields:
        try:
            field_name, alias = attr.split(':')
        except ValueError:
            field_name = alias = attr

        attr_names[field_name] = alias

    # Parse results of the homology search
    aln_totals = 0  #all hits to any database
    dis_totals = 0  #multiple annotations for single query
    hits = {}  #store matches and additional attributes
    for hits_file in args.hits:
        with open_io(hits_file) as in_h:
            if aln_format == "b6":
                reader = B6Reader(in_h)
                iter_obj = reader.iterate(header=specifiers)
            else:
                reader = HMMReader(in_h)
                iter_obj = reader.iterate()

            for entry in iter_obj:
                aln_totals += 1

                hit = Homolog()
                hit.parse_entry(entry, reverse_query=revq, domain=dom)

                query = hit.seqid

                if query in hits:
                    dis_totals += 1

                    if match_precedence == 'order':
                        # Preserve only best hit on first come basis
                        out_log("{}\t{}\t{}\n".format(hits[query].dbid, \
                                hit.dbid, 'order').encode('utf-8'))
                    else:
                        # Determine which match has the best alignment score
                        prev_score = hits[query].score

                        if prev_score >= hit.score:  #existing match is best
                            out_log("{}\t{}\t{}\n".format(hits[query].dbid, \
                                    hit.dbid, 'score').encode('utf-8'))
                        else:
                            out_log("{}\t{}\t{}\n".format(hit.dbid, \
                                    hits[query].dbid, 'score').encode('utf-8'))
                            hits[query] = hit
                else:
                    hits[query] = hit

    # Parse GFF file, adding annotations from homology search and relational db
    gff_totals = 0  #all features contained in GFF
    annot_totals = 0  #all features with annotation
    ftype_totals = 0  #features with type matching argument input
    no_map = 0  #track hits without an entry in the relational database
    no_id = 0  #track GFF3 entries without ID attribute to link to protein seq

    gff_reader = GFF3Reader(in_gff)
    for entry in gff_reader.iterate(parse_attr=True, headers=True):
        try:
            seq_id = entry.seqid
        except AttributeError:
            if entry.startswith('##'):
                write_io(out_h, "{}\n".format(entry))
                continue
            else:
                continue  #don't output comments

        gff_totals += 1

        try:
            feature_id = entry.attributes[id_attr]
        except KeyError:
            no_id += 1
            # Discard feature if unable to annotate and flag provided
            if not args.filter:
                write_io(out_h, entry.write())
            continue

        if type(feature_id) != string_type:
            line_number = gff_reader.current_line
            filename = os.path.basename(gff_reader.filename)
            print("error: {}, line {!s}: the identifier attribute '{}' cannot "
                  "have more than one value".format(filename, line_number, \
                  id_attr), file=sys.stderr)
            sys.exit(1)

        # Annotate features of a given type only
        if feature_type:
            if feature_type != entry.type:
                if not args.filter:
                    write_io(out_h, entry.write())
                continue
            else:
                ftype_totals += 1

        # clear existing attributes, if directed
        if args.clear_attrs:
            for attr in list(entry.attributes.keys()):
                if not attr[0].isupper():  #preserve reserved attributes
                    del(entry.attributes[attr])

        # Annotate features using attributes field
        try:
            # Query name should be same as value in identifier attribute
            hit = hits[feature_id]
        except KeyError:
            # Discard feature if unable to annotate and flag provided
            if args.filter:
                continue
        else:
            annot_totals += 1
            subject = hit.dbid

            entry.attributes['Alias'] = subject
            # Include additional information about a hit if a relational 
            # database provided
            if mapping:
                try:
                    sub_entry = mapping[subject]
                except KeyError:
                    no_map += 1
                else:
                    # Only include information from the requested fields
                    for field in all_fields:
                        attr_name = attr_names[field]

                        try:
                            entry_value = sub_entry[field]
                        except KeyError:
                            no_fields[field] += 1
                        else:
                            if entry_value:  #entry field might be blank
                                entry.attributes[attr_name] = entry_value

        # Add default to product attribute if it doesn't already exist
        if default_product:
            prod_name = attr_names['product']
            if prod_name not in entry.attributes:
                entry.attributes[prod_name] = default_product

        # Write annotated features to new GFF3 file
        write_io(out_h, entry.write())

    # Calculate and print statistics
    if no_map > 0:
        print("", file=sys.stderr)
        print("warning: there were {!s} alignments that did not have a "
              "corresponding mapping entry".format(no_map), file=sys.stderr)

    if no_id > 0:
        print("", file=sys.stderr)
        print("warning: there were {!s} features without an ID attribute tag "
              "linking them to a protein sequence".format(no_id), \
              file=sys.stderr)

    for map_field in no_fields:
        no_map_size = no_fields[map_field]
        if no_map_size > 0:
            print("", file=sys.stderr)
            print("warning: there were {!s} alignments that did not have an "
                  "entry with field '{}' in a mapping file"\
                  .format(no_map_size, map_field), file=sys.stderr)

    print("", file=sys.stderr)
    print("Alignments processed:", file=sys.stderr)
    print("  - alignment totals:\t{!s}".format(aln_totals), file=sys.stderr)
    print("  - conflicting annotations:\t{!s}"\
          .format(dis_totals), file=sys.stderr)
    print("Features processed:", file=sys.stderr)
    print("  - feature totals:\t{!s}".format(gff_totals), file=sys.stderr)
    if feature_type:
        print("  - of relevant type:\t{!s}".format(ftype_totals), \
              file=sys.stderr)
    print("  - with annotation:\t{!s}".format(annot_totals), \
          file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to annotate {!s} features"\
          .format(total_time, gff_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
