#! /usr/bin/env python
"""
Screen the results of a homology search using bit-score thresholds, 
alternative phenotype-conferring snps, or other scoring metrics.

Usage:
    screen_features [options] [-m in.json] [-s in.csv] [-o out.b6] in.b6

Required input is a tabular file of pairwise alignments (B6 format). Optional 
inputs are a tabular SNP file containing information on alternative phenotypes 
and a relational database in JSON format containing an appropriate scoring 
threshold (e.g. bitscore thresholds).

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension 
to the output file name (e.g. .gz, .bz2). Use /dev/stdin as the file name to 
indicate that input is from standard input (stdin). Similarly, leave off 
'--out' to direct output to standard output (stdout).

Copyright:

    screen_features Screens homology search results using various metrics
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
from bio_utils.iterators import b6_iter
import json
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.1.1"


def parse_snps_file(in_h, line=None):
    """
    Parse an SNPs file. Each line should have six columns, corresponding to 
    accession number, model type, parameter type, position of SNP, wildtype 
    base composition, and mutant base composition. The header will be 
    ignored, if present.
    """

    append = list.append
    join = str.join
    strip = str.strip

    if line is None:
        line = next(in_h)  # Read header

    # Check if input is text or bytestream
    if (isinstance(line, bytes)):
        def next_line(i):
            return next(i).decode('utf-8')

        line = line.decode('utf-8')
    else:
        next_line = next

    # Check if first line is header
    if line.startswith('#'):
        split_line = strip(next_line(in_h)).split('\t')
    else:
        split_line = strip(line).split('\t')

    snp_dict = {}
    try:

        while True:  #loop until StopIteration raised

            try:
                acc, mt, pt, position, wildtype, mutant = split_line
            except ValueError:
                raise FormatError("a tabular SNP file is required with six "
                                  "columns, corresponding to accession #, "
                                  "model type, parameter type, position, and "
                                  "wildtype and mutatant alleles")

            try:
                snp_dict[acc].append((position, wildtype, mutant))
            except KeyError:
                snp_dict[acc] = [(position, wildtype, mutant)]

            split_line = strip(next_line(in_h)).split('\t')

    except StopIteration:

        return snp_dict


def screen_aln_quality(hit, e=10, ident=0, cov=0, score=0):
    """
    Screen reference match based on alignment quality metrics.
    """
    condition = (hit.evalue <= e and hit.perc_identical >= ident and \
                 hit.align_len >= cov and hit.bit_score >= float(score))

    if condition:
        return True
    else:
        return False


def screen_snp(hit, acc, snps, header, only_snp=False):
    """
    Screen reference match for alternative phenotype-conferring SNPs using 
    a dictionary of accession number, SNP value pairs.
    """

    try:
        mutations = snps[acc]
    except KeyError:
        print("warning: '{}' not found in the database of SNPs".format(acc), \
              file=sys.stderr)

        if only_snp:
            return False
        else:
            return True

    else:

        try:
            qseq_index = header.index('qseq') - 12  #minus required indices
            qseq = hit.add_specs[qseq_index]
        except KeyError:
            print("error: the format specifier 'qseq' is required for "
                  "secondary screening of SNPs".format(), file=sys.stderr)
            sys.exit(1)

        try:
            sseq_index = header.index('sseq') - 12
            sseq = hit.add_specs[sseq_index]
        except KeyError:
            print("error: the format specifier 'sseq' is required for "
                  "secondary screening of SNPs".format(), file=sys.stderr)
            sys.exit(1)

        start = int(hit.subject_start)
        end = int(hit.subject_end)

        for mutation in mutations:
            pos, wt, sub = mutation
            pos = int(pos)

            if pos not in list(range(start, end + 1)):
                return False

            else:
                rel_pos = pos - start  #position of snp relative to alignment

                # Handle any gaps
                i = 0
                while i != (rel_pos + 1):
                    if sseq[i] == '-':
                        rel_pos += 1
                        i += 1
                    else:
                        i += 1

                # Handle database inconsistencies
                if sseq[rel_pos] != wt:
                    print("warning: possible inconsistency in database: "
                          "subject {} with residue {} at SNP position {!s} "
                          "does not match wildtype residue {}"\
                          .format(hit.subject, sseq[rel_pos], pos, wt), \
                          file=sys.stderr)
                    return False

                if qseq[rel_pos] == sub:
                    return True
                else:
                    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('b6',
        metavar='in.b6',
        action=Open,
        mode='rb',
        default=sys.stdin,
        help="input best-hit alignments in B6/M8 format. The alignment format "
             "should be provided to -s/--specifiers if other than the default "
             "BLAST+ tabular format")
    parser.add_argument('-o', '--out',
        metavar='out.b6',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output screened alignments in B6/M8 format [default: output to "
              "stdout]")
    parser.add_argument('-s', '--specifiers',
        dest="format",
        action=ParseSeparator,
        sep=",",
        help="input ordered field specifiers [default: qaccver, saccver, "
             "pident, length, mismatch, gapopen, qstart, qend, sstart, send, "
             "evalue, bitscore]. The additional field specifiers 'qseq' and "
             "'sseq' are required when screening for SNPs")
    screen_method = parser.add_argument_group(title="screening parameters")
    screen_method.add_argument('-e', '--evalue',
        type=float,
        default=10,
        help="maximum expect value allowed for a hit to be retained [default: "
             "10]")
    screen_method.add_argument('-m', '--mapping',
        metavar='mapping.json',
        dest='map_file',
        action=Open,
        mode='rb',
        help="input relational database in JSON format containing "
             "model-specific bit-score or alternative scoring thresholds")
    screen_method.add_argument('-f', '--score-field',
        dest="score_field",
        help="field in the relational database to extract the alignment "
             "quality thresholds from")
    screen_method.add_argument('-c', '--category-field',
        dest="cat_field",
        default="gene_family",
        help="field in the relational database to extract the category "
             "accessions from [default: gene_family]")
    screen_method.add_argument('-i', '--identity',
        type=float,
        default=0,
        help="minimum percent identity required to retain a hit [default: 0]")
    screen_method.add_argument('-l', '--length',
        dest='aln_len',
        type=int,
        default=0,
        help="minimum alignment length required to retain a hit [default: 0]. "
             "Should ideally be used with --identity, as alignment length is "
             "not useful on its own")
    screen_method.add_argument('-a', '--alleles',
        metavar='snps.csv',
        dest="snps",
        action=Open,
        mode='rb',
        help="tabular file containing SNPs (position, along with the wildtype "
             "and mutant alleles) that confer alternative phenotypes to an "
             "organism. This step can only be used with -m/--mapping and will "
             "be performed last if any additional screening for alignment "
             "quality is also specified")
    screen_method.add_argument('--snp',
        dest='only_snp',
        action='store_true',
        help="only output matches if found in the database of SNPs [default: "
             "False]")
    parser.add_argument('--defaults',
        dest='default_format',
        action='store_true',
        help="only output default BLAST+ specifiers [default: output all]")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if not (args.map_file or args.evalue or args.snps or args.aln_len or \
            args.identity):
        parser.error("error: one or more of the following arguments are "
                     "required: -s/--snps, -e/--evalue, -m/--mapping, "
                     "-i/--identity, or -l/--length")

    if args.snps and not args.map_file:
        parser.error("error: -m/--mapping required when -s/--snps is supplied")

    if args.score_field and not args.map_file:
        parser.error("error: -m/--mapping required when -f/--score-field is "
                     "supplied")


    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('screen_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)


    # Track program run-time
    start_time = time()


    # Assign variables based on user inputs
    out_h = args.out.write

    if args.map_file:
        mapping = json.load(args.map_file)

    snps = parse_snps_file(args.snps) if args.snps else None

    only_snp = args.only_snp
    specifiers = args.format
    score_field = args.score_field
    cat_field = args.cat_field
    e_thresh = args.evalue
    id_thresh = args.identity
    len_thresh = args.aln_len
    default_only = True if args.default_format else False


    # Screen hits for alignment quality and/or mutant alleles
    passed_total = 0
    for totals, hit in enumerate(b6_iter(args.b6, header=specifiers)):

        subject = hit.subject

        if args.map_file:
            try:
                sub_entry = mapping[subject]
            except KeyError:
                print("error: subject id '{}' not found in the relational "
                      "database".format(subject), file=sys.stderr)
                sys.exit(1)

            if score_field:
                try:
                    score = sub_entry[score_field]
                except KeyError:
                    print("error: field '{}' not found in the relational "
                          "database".format(score_field), file=sys.stderr)
                    sys.exit(1)

            else:
                score = 0

        else:
            score = 0

        q_pass = screen_aln_quality(hit, e=e_thresh, score=score, \
                                    ident=id_thresh, cov=len_thresh)

        if q_pass:
            if snps:
                try:
                    cat_acc = sub_entry[cat_field]
                except UnboundLocalError:
                    print("error: no entry in the relational database for {}"\
                          .format(subject), file=sys.stderr)
                    sys.exit(1)

                s_pass = screen_snp(hit, cat_acc, snps, header=specifiers, \
                                    only_snp=only_snp)

                if s_pass:
                    passed_total += 1
                    out_h(hit.write(defaults=default_only))

            else:
                passed_total += 1
                out_h(hit.write(defaults=default_only))


    # Output screening statistics
    print("Total number of matches:\t{!s}".format(totals), \
          file=sys.stderr)
    print("  - matches that passed screening:\t{!s}\n".format(passed_total), \
          file=sys.stderr)
    

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to screen {!s} matches\n"\
          .format(total_time, totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
