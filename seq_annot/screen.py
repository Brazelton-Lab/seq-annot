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

from arandomness.argparse import Open, ParseSeparator
import argparse
from bio_utils.iterators import B6Reader
import json
import os
import re
from seq_annot.db import load_dbs
from seq_annot.seqio import write_io
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.4.5"


def screen_aln_quality(hit, evalue=10, identity=0, length=0, score=0):
    """
    Screen results of a homology search based on alignment quality metrics.

    Args:
        hit (B6Entry): Class containing all B6/M8 data

        evalue (float): expect value maximum [default: 10]

        identity (float): minimum percentage identity shared between query and 
            subject [default: 0]

        length (int): minimum length of the alignment [default: 0]

        score (int): bitscore threshold [default: 0]

    Returns:
        bool: True if hit passes screening else False
    """
    passed = (hit.evalue <= evalue and hit.identity >= identity and \
                 hit.length >= length and hit.bitscore >= float(score))

    if passed:
        return True
    else:
        return False


def screen_snp(hit, snps, only_snp=False):
    """
    Screen results of a homology search for alternative phenotype-conferring 
    SNPs.

    Args:
        hit (B6Entry): Class containing all B6/M8 data 

        snps (list): List of SNPs with format <wild-type residue>
            <position><mutant residue> (e.g. A256G)

        only_snps (bool): Returns False when SNP information is not available
            for the query subject

    Returns:
        bool: True if hit passes screening else False
    """
    SNP = re.compile("([A-Z])([0-9]+)([A-Z])")

    # Handle case where no SNPs are available for the given subject
    if not snps:
        if only_snp:
            return False
        else:
            return True

    try:
        qseq = hit.custom_fs['qseq']
    except KeyError:
        print("error: the format specifier 'qseq' is required for "
              "secondary screening of SNPs".format(), file=sys.stderr)
        sys.exit(1)

    try:
        sseq = hit.custom_fs['sseq']
    except KeyError:
        print("error: the format specifier 'sseq' is required for "
              "secondary screening of SNPs".format(), file=sys.stderr)
        sys.exit(1)

    start = int(hit.subject_start)
    end = int(hit.subject_end)

    match = False
    for snp in snps:
        try:
            wt, pos, sub = SNP.search(snp).groups()
        except AttributeError:
            print("error: SNP {} from subject {} is formatted incorrectly. "
                  "Please see the help message for formatting requirements."\
                  .format(snp, hit.subject))
            sys.exit(1)

        pos = int(pos)

        if pos not in list(range(start, end + 1)):
            continue

        else:
            rel_pos = pos - start  #position of snp relative to alignment

            # Handle gaps
            i = 0
            while i != (rel_pos + 1):
                if sseq[i] == '-':
                    rel_pos += 1
                    i += 1
                else:
                    i += 1

            if qseq[rel_pos] == sub:
                match = True
                    
    return match


def do_nothing(*args):
    pass


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
        mode='wb',
        default=sys.stdout,
        help="output screened alignments in B6/M8 format [default: output to "
              "stdout]")
    parser.add_argument('-d', '--discards',
        dest='discards',
        metavar='out.csv',
        action=Open,
        mode='wb',
        help="output discards to file")
    parser.add_argument('-s', '--specifiers',
        dest="format",
        action=ParseSeparator,
        sep=",",
        default=["qaccver", "saccver", "pident", "length", "mismatch",
                 "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                 "bitscore"],
        help="input ordered field specifiers [default: qaccver, saccver, "
             "pident, length, mismatch, gapopen, qstart, qend, sstart, send, "
             "evalue, bitscore]. The additional field specifiers 'qseq' and "
             "'sseq' are required when screening for SNPs")
    parser.add_argument('-m', '--mapping',
        metavar='mapping.json',
        dest='map_file',
        action=Open,
        mode='rb',
        help="input relational database in JSON format. Required when "
             "screening for SNPs or alignment quality using bitscore or "
             "model-specific scoring thresholds")
    screen_method = parser.add_argument_group(title="screening parameters")
    screen_method.add_argument('-e', '--evalue',
        metavar='VALUE',
        type=float,
        default=10,
        help="maximum expect value allowed for a hit to be retained [default: "
             "10]")
    screen_method.add_argument('-b', '--bitscore',
        metavar='FIELD',
        dest="score_field",
        help="field in the relational database containing scoring thresholds "
             "that will be used to screen alignment quality. Argument must be "
             "used in conjunction with -m/--mapping")
    screen_method.add_argument('-i', '--identity',
        metavar='PERCENT',
        type=float,
        default=0,
        help="minimum percent identity required to retain a hit [default: 0]")
    screen_method.add_argument('-l', '--length',
        metavar='LENGTH',
        dest='aln_len',
        type=int,
        default=0,
        help="minimum alignment length required to retain a hit [default: 0]")
    screen_method.add_argument('-a', '--alleles',
        metavar='FIELD',
        dest="snp_field",
        help="field in the relational database containing SNPs conferring an "
             "alternative phenotype to an organism. The field should contain "
             "a list of one or more SNPS with format <wild-type "
             "residue><position><mutant residue> (e.g. A254G). Argument must "
             "be used in conjunction with -m/--mapping. Screening for SNPs "
             "will be performed after alignment quality screening")
    output_control = parser.add_argument_group(title="output control options")
    output_control.add_argument('--snp',
        dest='only_snp',
        action='store_true',
        help="only output matches if the matches have corresponding SNPs "
             "[default: False]")
    output_control.add_argument('--defaults',
        dest='default_format',
        action='store_true',
        help="only output default BLAST+ specifiers [default: output all]")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if not (args.evalue or args.snp_field or args.aln_len or args.identity or \
            args.score_field):
        parser.error("error: one or more of the following arguments is "
                     "required: -a/--alleles, -e/--evalue, -b/--bitscore, "
                     "-i/--identity, or -l/--length")

    if args.snp_field and not args.map_file:
        parser.error("error: -m/--mapping required when -a/--alleles is "
                     "supplied")

    if args.score_field and not args.map_file:
        parser.error("error: -m/--mapping required when -b/--bitscore is "
                     "supplied")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('screen_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    out_h = args.out

    if args.discards:
        out_d = args.discards.write
        out_d("#Discarded\tReason\n".encode('utf-8'))
    else:
        out_d = do_nothing

    mapping = load_dbs(args.map_file) if args.map_file else None

    specifiers = args.format
    score_field = args.score_field
    snp_field = args.snp_field
    e_thresh = args.evalue
    id_thresh = args.identity
    len_thresh = args.aln_len
    default_only = True if args.default_format else False
    only_snp = args.only_snp

    # Screen hits for alignment quality and/or mutant alleles
    passed_total = 0  #hits passed screening
    failed_qual = 0  #hits failed due to alignment quality
    failed_snp = 0  #hits failed due to SNP screening
    failed_bh = 0  #hits failed due to best-hit criteria
    aln_totals = 0  #all alignments in B6 file

    prev_hit = ''
    b6_reader = B6Reader(args.b6)
    for hit in b6_reader.iterate(header=specifiers):
        aln_totals += 1

        subject = hit.subject

        if mapping:
            try:
                sub_entry = mapping[subject]
            except KeyError:
                line_number = b6_reader.current_line
                filename = os.path.basename(b6_reader.filename)
                print("error: {}, line {!s}: subject id {} not found in the "
                      "database".format(filename, line_number, subject), 
                      file=sys.stderr)
                sys.exit(1)

        if score_field:
            try:
                score = sub_entry[score_field]
            except KeyError:
                line_number = b6_reader.current_line
                filename = os.path.basename(b6_reader.filename)
                print("error: {}, line {!s}: {} has no field named {} in "
                      "database".format(filename, line_number, subject, 
                      score_field), file=sys.stderr)
                sys.exit(1)
            else:
                try:
                    score = float(score)
                except ValueError:  #no value in field
                    print("", file=sys.stderr)
                    print("warning: no value for {} in field {}. Setting "
                          "scoring threshold to 0\n".format(score_field, \
                          subject, hit.query), file=sys.stderr)
                    score = 0

        else:
            score = 0

        q_pass = screen_aln_quality(hit, evalue=e_thresh, score=score, \
                                    identity=id_thresh, length=len_thresh)

        if q_pass:
            if snp_field:
                try:
                    snps = sub_entry[snp_field]
                except KeyError:
                    line_number = b6_reader.current_line
                    filename = os.path.basename(b6_reader.filename)
                    print("error: {}, line {!s}: {} has no field named {} in "
                          "database".format(filename, line_number, subject, 
                          snp_field), file=sys.stderr)
                    sys.exit(1)

                s_pass = screen_snp(hit, snps, only_snp=only_snp)

                if s_pass:
                    passed_total += 1
                    write_io(out_h, hit.write(default=default_only))
                else:
                    failed_snp += 1
                    out_d("{}\tSNP screen\n".format(hit.query).encode('utf-8'))

            else:
                passed_total += 1
                write_io(out_h, hit.write(default=default_only))

        else:
            failed_qual += 1
            out_d("{}\talignment quality\n".format(hit.query).encode('utf-8'))

    # Output screening statistics
    print("", file=sys.stderr)
    print("Alignments processed:", file=sys.stderr)
    print("  - alignment totals:\t{!s}".format(aln_totals), file=sys.stderr)
    print("  - passed screening:", file=sys.stderr)
    print("    - passed totals:\t{!s}".format(passed_total), file=sys.stderr)
    print("  - failed screening:", file=sys.stderr)
    print("    - failed totals:\t{!s}".format(aln_totals - \
          passed_total), file=sys.stderr)
    if e_thresh or len_thresh or id_thresh or score_field:
        print("    - from alignment quality:\t{!s}".format(failed_qual), \
              file=sys.stderr)
    if snp_field:
        print("    - from secondary SNP screening:\t{!s}".format(failed_snp), \
              file=sys.stderr)
 
    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("", file=sys.stderr)
    print("It took {:.2e} minutes to screen {!s} alignments"\
          .format(total_time, aln_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
