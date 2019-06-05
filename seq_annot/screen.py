#! /usr/bin/env python
"""
Screen the results of a homology search using alignment scoring metrics such 
as bit-score thresholds, alternative phenotype-conferring snps, sequence 
similarity percentage, etc.

Required input is a pairwise alignments file in B6 format. Optional inputs 
are a tabular SNP file containing information on alternative phenotypes and 
a relational database in JSON format containing an appropriate scoring 
threshold (e.g. bitscore thresholds).

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension 
to the output file name (e.g. .gz, .bz2). Use "-" indicate that input is from 
standard input (stdin). Similarly, leave off "--out" to direct output to 
standard output (stdout).

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

import argparse
from bio_utils.iterators import B6Reader
import json
import os
import re
from seq_annot.argparse import *
from seq_annot.db import load_dbs
from seq_annot.seqio import *
import sys
import textwrap
from time import time

SNP = re.compile("([A-Z])([0-9]+)([A-Z])")

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.5.3"


class SubjectHSPs:
    """
    Attributes:
        query (str): query ID (sequence aligned with)

        subject (str): subject ID (sequence aligned to)

        hsp (list): list containing B6Entry objects for each high-scoring 
            segment pairs of a subject sequence
    """

    def __init__(self):
        """Initialize variables to store B6/M8 entry data"""

        self.query = None
        self.subject = None
        self.hsp = []

    def screen(self, e_thresh=10, i_thresh=0, sim_thresh=0, cov_thresh=0, \
        len_thresh=0, score_thresh=0):
        """Screen results of a homology search based on alignment quality 
        metrics.

        Args:
            e_thresh (float): expect value maximum [default: 10]

            i_thresh (float): minimum percentage identity shared between query 
                and subject [default: 0]

            sim_thresh (float): minimum percentage sequence similarity shared 
                between query and subject [default: 0]

            cov_thresh (float): minimum percentage coverage shared between query 
                and subject [default: 0]

            len_thresh (int): minimum length of the alignment [default: 0]

            score_thresh (int): bitscore threshold [default: 0]

        Returns:
            bool: True if subject passes screening else False
        """
        e = []
        i = []
        l = []
        b = []
        for hit in self.hsp:
            e.append(hit.evalue)
            i.append(hit.identity)
            l.append(hit.length)
            b.append(hit.bitscore)

        evalue = min(e)
        bitscore = max(b)
        length = max(l)
        identity = max(i)

        if cov_thresh:
            try:
                cov = int(self.hsp[0].custom_fs['qcovs'])
            except KeyError:
                raise FormatError("field specifier 'qcovs' required when "
                    "coverage threshold provided")
            except IndexError:
                raise AlignmentNotFound("")
            except ValueError:
                raise FormatError("'qcovs' must be an integer value")
        else:
            cov = 1

        if sim_thresh:
            sim = sum([i[pos] * l[pos] for pos in range(len(i))]) / sum(l)
        else:
            sim = 1

        try:
            passed = (evalue <= e_thresh and identity >= i_thresh and \
                length >= len_thresh and bitscore >= float(score_thresh) and \
                cov >= cov_thresh and sim >= sim_thresh)
        except TypeError:
            raise FormatError("minimum format specifiers required for "
                "screening are 'evalue, pident, length, and bitscore'")

        return(passed)

    def snp(self, snps):
        """Screen results of a homology search for alternative 
        phenotype-conferring SNPs.

        Args:
            snps (list): List of SNPs with format <wild-type residue>
                <position><mutant residue> (e.g. A256G)

        Returns:
            bool: True if hit passes screening else False
        """
        # Handle case where no SNPs are available for the given subject
        if not snps:
            return(True)

        match = False
        for hit in self.hsp:
            start = int(hit.subject_start)
            end = int(hit.subject_end)

            try:
                qseq = hit.custom_fs['qseq']
            except KeyError:
                raise FormatError("format specifier 'qseq' is required for "
                    "secondary screening of SNPs", file=sys.stderr)

            try:
                sseq = hit.custom_fs['sseq']
            except KeyError:
                raise FormatError("format specifier 'sseq' is required for "
                    "secondary screening of SNPs", file=sys.stderr)

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
                        return(True)
                    
        return match

    def write(self, handle, defaults=False):
        """
        Args:
            handle (File): File handle of output file

            defaults (bool): True if output format should be in standard B6 
                format, otherwise output as is
        """
        for hit in self.hsp:
            write_io(handle, hit.write(default=defaults))


def screen_query(hit, out_h, out_d, metrics, defaults=False):
    """Run quality and SNP screening on hits
    """
    passed = 0
    failed_qual = 0
    failed_snp = 0

    if not hit.query:
        return(passed, failed_qual, failed_snp)

    quals = metrics[0:-1]
    snps = metrics[-1]

    q_pass = hit.screen(*quals)

    if q_pass:
        s_pass = hit.snp(snps)

        if s_pass:
            passed = 1
            hit.write(out_h, defaults)
        else:
            failed_snp = 1
            out_d("{}\tSNP screen\n".format(hit.query).encode('utf-8'))
    else:
        failed_qual = 1
        out_d("{}\talignment quality\n".format(hit.query).encode('utf-8'))

    return(passed, failed_qual, failed_snp)

def restricted_float(x):
    """Raise error if argument value out of restricted range
    """
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("{!s} not in range [0.0, 1.0]"\
            .format(x))
    return x

def do_nothing(*args):
    pass

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('b6',
        metavar='in.b6',
        help="input best-hit alignments in B6/M8 format. The alignment format "
            "should be provided to -s/--specifiers if other than the default "
            "BLAST+ tabular format. Use '-' to accept input from stdin")
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
    parser.add_argument('-f', '--specifiers',
        dest="format",
        action=ParseSeparator,
        sep=",",
        default=["qaccver", "saccver", "pident", "length", "mismatch",
            "gapopen", "qstart", "qend", "sstart", "send", "evalue", 
            "bitscore"],
        help="input ordered field specifiers [default: qaccver, saccver, "
            "pident, length, mismatch, gapopen, qstart, qend, sstart, send, "
            "evalue, bitscore]. The additional field specifiers 'qseq' and "
            "'sseq' are required when screening for SNPs and the specifier "
            "'qcovs' is required when screening for coverage")
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
        metavar='[0.0-1.0]',
        type=restricted_float,
        default=0,
        help="minimum identity threshold (0-1) required to retain a hit "
            "[default: 0]")
    screen_method.add_argument('-s', '--similarity',
        metavar='[0.0-1.0]',
        dest="similarity",
        type=restricted_float,
        default=0,
        help="minimum alignment similarity threshold (0-1) required to retain "
            "a hit. Sequence similarity is calculated as the sum of identical "
            "matches divided by the sum of alignment lengths for all HSPs per "
            "subject [default: 0]")
    screen_method.add_argument('-c', '--coverage',
        metavar='[0.0-1.0]',
        dest="coverage",
        type=restricted_float,
        default=0,
        help="minimum coverage per subject threshold (0-1) required to retain "
            "a hit [default: 0]")
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
    output_control.add_argument('--defaults',
        dest='default_format',
        action='store_true',
        help="only output default BLAST+ specifiers [default: output all]")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if args.snp_field and not args.map_file:
        parser.error("error: -m/--mapping required when -a/--alleles is "
                     "supplied")

    if args.score_field and not args.map_file:
        parser.error("error: -m/--mapping required when -b/--bitscore is "
                     "supplied")
    
    if args.coverage and not ("qcovs" in args.format):
        parser.error("error: field specifier 'qcovs' required when argument "
            "-c/--coverage is supplied")

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

    if args.b6 == '-':
        b6 = sys.stdin
    else:
        b6 = open_io(args.b6, 'rb')

    mapping = load_dbs(args.map_file) if args.map_file else None

    specifiers = args.format
    score_field = args.score_field
    snp_field = args.snp_field
    e_thresh = args.evalue
    id_thresh = args.identity
    len_thresh = args.aln_len
    sim_thresh = args.similarity * 100
    cov_thresh = args.coverage * 100
    default_only = True if args.default_format else False

    # Screen hits for alignment quality and/or mutant alleles
    passed_total = 0  #hits passed screening
    failed_quals = 0  #hits failed due to alignment quality
    failed_snps = 0  #hits failed due to SNP screening
    aln_totals = 0  #all alignments in B6 file
    hit_totals = 0  #all hits to a subject sequence in B6 file

    prev_hit = ''
    b6_reader = B6Reader(b6)
    hit = SubjectHSPs()
    metrics = []
    for entry in b6_reader.iterate(header=specifiers):
        aln_totals += 1

        subject = entry.subject
        query = entry.query
        if query != hit.query:  #new query encountered
            hit_totals += 1

            # Process previous query
            totals = screen_query(hit, out_h, out_d, metrics, default_only)
            passed_total += totals[0]
            failed_quals += totals[1]
            failed_snps += totals[2]

            # Reset variables for new query
            hit = SubjectHSPs()
            hit.query = query
            hit.subject = subject

            if mapping:
                try:
                    sub_entry = mapping[subject]
                except KeyError:
                    line_number = b6_reader.current_line
                    filename = os.path.basename(b6_reader.filename)
                    raise FormatError("{}, line {!s}: subject id {} not found "
                        "in database".format(filename, line_number, subject), 
                        file=sys.stderr)

            if score_field:
                try:
                    score = sub_entry[score_field]
                except KeyError:
                    line_number = b6_reader.current_line
                    filename = os.path.basename(b6_reader.filename)
                    raise FormatError("{}, line {!s}: {} has no field named {} "
                        "in database".format(filename, line_number, subject, 
                        score_field), file=sys.stderr)
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

            if snp_field:
                try:
                    snps = sub_entry[snp_field]
                except KeyError:
                    line_number = b6_reader.current_line
                    filename = os.path.basename(b6_reader.filename)
                    raise FormatError("{}, line {!s}: {} has no field named {} "
                        "in database".format(filename, line_number, subject, 
                        snp_field), file=sys.stderr)
            else:
                snps = []

            metrics = [e_thresh, id_thresh, sim_thresh, cov_thresh, \
                len_thresh, score, snps]

        # Add HSP to hit
        hit.hsp.append(entry)

    # Process final query
    totals = screen_query(hit, out_h, out_d, metrics, default_only)
    passed_total += totals[0]
    failed_quals += totals[1]
    failed_snps += totals[2]

    # Output screening statistics
    print("", file=sys.stderr)
    print("Alignments processed:", file=sys.stderr)
    print("  - hit totals:\t{!s}".format(hit_totals), file=sys.stderr)
    print("  - alignment totals:\t{!s}".format(aln_totals), file=sys.stderr)
    print("Screening results:", file=sys.stderr)
    print(" - passed totals:\t{!s}".format(passed_total), file=sys.stderr)
    print(" - failed totals:\t{!s}".format(failed_quals + failed_snps), \
        file=sys.stderr)
    if e_thresh or len_thresh or id_thresh or score_field or cov_thresh \
        or sim_thresh:
        print("    - from alignment quality:\t{!s}".format(failed_quals), \
              file=sys.stderr)
    if snp_field:
        print("    - from secondary SNP screening:\t{!s}".format(failed_snps),\
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
