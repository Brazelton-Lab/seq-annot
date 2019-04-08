#!/usr/bin/env python
"""
Estimate abundance for genomic features.

Required inputs are a GFF3 file of annotated features and an alignments file 
in SAM or BAM format. Optional input is one or more relational databases 
containing features mapped to feature categories, such as genes to gene
families. Abundances will be calculated for feature categories in place of 
features if the appropriate arguments are supplied. Output is a tabular file 
of estimated feature abundances.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, provide the appropriate flag based on 
the desired compression algorithm. Standard input (stdin) can be redirected to 
one of the positional arguments if supplied with '-'.

Copyright:

    count_features Counts coverage of GFF entries from a SAM/BAM file
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
from __future__ import division

from arandomness.argparse import Open,ParseSeparator
import argparse
import HTSeq
import itertools
import json
import os
from seq_annot.seqio import open_io, write_io
from statistics import mean
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Beta"
__version__ = '1.5.6'


class UnknownChrom(Exception):
   pass


def invert_strand(iv):
    """Return inverted strand
    """
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        print("error: unknown strand value {}".format(iv2), file=sys.stderr)
        sys.exit(1)
    return iv2


def compute_effective_length(length, fld=None, mu=600):
    """Compute the effective length of a feature. The effective length is 
    defined as the number of possible start positions for a fragment within 
    a feature, given that the fragment must be contained entirely within its 
    boundaries.

    formula:
        effLen_i = featureLen_i - mu_i + 1

    featureLen - length of the feature
    mu - mean of the FDL for all fragments of length < featureLen
    FLD - fragment length (insert-size) distribution
    """
    mean_fld = mean([i for i in fld if i < length]) if fld else mu
    return length - mean_fld + 1


def scale_abundance_none(counts, *args):
    """Return untransformed coverage values.
    """
    return counts


def scale_abundance_fpk(counts, length, *args):
    """transform feature counts to fragements per kilobase of feature.

    formula:
        fpk_i = counts_i / ( effLen_i/1000 )

    counts - number of reads mapped to a feature (i.e. a gene sequence)
    effLen - effective length of a feature
    """
    return counts / (length / 1000)


def scale_abundance_fpkm(counts, length, sfactor, *args):
    """Transform feature counts to fragements per kilobase of feature per 
    million mapped fragements.

    formula:
        fpkm_i = counts_i / ( effLength_i/1,000 * totalNumReads_i/1,000,000 )

    counts - number of reads mapped to a feature (i.e. a gene sequence)
    effLen - effective length of a feature
    totalNumReads - total number of mapped reads per sample

    reference DOI: doi:10.1038/nmeth.1226
    """
    return counts / ((length / 1000) * (sfactor / 1000000))


def scale_abundance_tpm(counts, length, sfactor, *args):
    """Transform feature counts to transcripts per million.

    formula:
       tpm_i = countPerBase_i * ( 1/totalCountsPerBase ) * 1,000,000

    countsPerBase - number of counts per base ( counts/effLen )
    counts - number of reads mapped to a feature (i.e. a gene sequence)
    effLen - effective length of a feature
    totalCountsPerBase - sum of all counts per base rates

    reference DOI: doi:10.1007/s12064-012-0162-3
    """
    return (counts / length) * (1 / sfactor) * 1000000


def scale_abundance_prop(counts, length, sfactor, *args):
    """Transform feature counts using custom method.

    formula:
        prop_i = countPerBase_i / totalCountsPerBase * scaleFactor

    countsPerBase - number of counts per base ( counts/effLen )
    counts - number of reads mapped to a feature (i.e. a gene sequence)
    effLen - effective length of a feature
    totalCountsPerBase - sum of all counts per base rates
    """
    return (counts / length) * sfactor


def return_first_arg(first, *args):
    return first


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('alignment_file',
        metavar='in.aln',
        help="input alignment file in SAM or BAM format. Use '-' to indicate "
             "that input should be taken from standard input (stdin)")
    parser.add_argument('feature_file',
        metavar='in.gff3',
        help="input feature annotation file in GFF3 format. Use '-' to indicate "
             "that input should be taken from standard input (stdin)")
    parser.add_argument('-m', '--mapping', 
        metavar='in.json',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input one or more relational databases, in JSON format, "
             "containing features mapped to feature categories, such as genes "
             "to gene families or exons to genes. Abundance estimates for the "
             "given feature category will be reported in place of features. "
             "Multiple input files can be provided by separating them with a "
             "comma and no spaces")
    parser.add_argument('-o', '--outpref',
        type=str,
        metavar='PREFIX',
        dest='outpref',
        default='sample',
        help="prefix for the output tabular files containing feature abundance "
             "estimates [default: sample]. File names will be appended with "
             "the units, file format, and compression algorithm, if relevant "
             "[e.g. sample.counts.csv.gz]")
    parser.add_argument('-f', '--format', 
        metavar='FORMAT', 
        dest='aformat',
        choices=['bam', 'sam'], 
        default='bam',
        help="input alignment file format [default: bam]. Options are 'sam' "
             "or 'bam'")
    parser.add_argument('-c', '--category',
        metavar='FIELD',
        dest='category',
        help="field in the relational database representing how features "
             "are categorized")
    parser.add_argument('-s', '--sorting', 
        metavar='ORDER',
        dest='order',
        choices=["position", "name"], 
        default='position',
        help="alignment file sorting scheme. Options are 'position' and "
             "'name' [default: position]. Alignments must be pre-sorted "
             "either by position/coordinates or by read name. This option "
             "will be ignored for single-end reads")
    parser.add_argument('-t', '--type', 
        metavar='TYPE', 
        dest='ftype',
        default='CDS',
        help="feature type (3rd column in GFF file) to estimate abundance for "
             "[default: CDS]. All features of other type will be ignored")
    parser.add_argument('-a', '--attr',
        metavar='ATTRIBUTE',
        default="Name",
        help="GFF attribute to use as the ID for the calculated abundances "
             "[default: 'Name']. This value will also be used as the search "
             "ID in the relational database, if provided")
    parser.add_argument('-e', '--mode', 
        metavar='MODE',
        choices=["union", "intersection-strict", "intersection-nonempty"],
        default="union",
        help="mode for handling different alignment scenarios. Options are "
             "'union', 'intersection-strict', and 'intersection-nonempty' "
             "[default: union]. The modes will count alignments differently "
             "depending on whether a read/pair overlaps more than one feature "
             "or only partially aligns to a single feature. The most "
             "inclusive mode is 'union' when given with the nonunique flag, "
             "and the least inclusive is 'intersection-strict'")
    parser.add_argument('-u', '--units', 
        metavar='UNITS',
        dest='norm',
        action=ParseSeparator,
        sep=',',
        default='counts',
        help="comma-separated list of units to output abundance estimates in "
             "[default: counts]. Options are 'counts', 'fpk' (fragments per "
             "kilobase of feature), 'fpkm' (fragements per kilobase of "
             "feature per million mapped fragments), 'tpm' "
             "(transcripts/fragments per million), 'prop', and 'custom'. If "
             "other than 'counts', features will be normalized by recruitment "
             "length, which will be calculated from the start and end fields "
             "of the GFF3 file. This is the sole normalization method used "
             "when transforming counts to FPK, and is useful to correct for "
             "differences in feature lengths within a sample. In addition to "
             "feature length, FPKM and TPM attempt to account for differences "
             "between samples in sequencing effort. An advantage of TMP over "
             "FPKM is that TPM is a proportional measurement, making it "
             "easier to identify the extent that the relative 'importance' of "
             "a given feature changes between samples. A custom transformation "
             "can also be performed when used with the -k/--coeff argument, in "
             "which case the length normalized proportion of a feature will be "
             "multiplied by the provided scaling factor.")
    parser.add_argument('-k', '--coeff',
        metavar='MUL',
        dest='sfactor',
        type=float,
        default=1,
        help="multiplier to use when 'custom' is given to -u/--units "
             "[default: 1]")
    parser.add_argument('-q', '--qual', 
        metavar='THRESH',
        dest='minqual',
        type=int, 
        default=2,
        help="skip all reads with alignment quality lower than the threshold "
             "[default: 2]")
    parser.add_argument('-b', '--buffer', 
        metavar='BYTES',
        dest='buffer_size',
        type=int,
        default=3145728,
        help="buffer size for paired reads in the alignment file if sorted by "
             "position [default: 3145728 (3GB)]. This value should be "
             "increased if memory issues are encountered")
    parser.add_argument('--cdna',
        dest='transcripts',
        action='store_true',
        help="sequences represent cDNA [default: False]. Whether sequences are "
             "from gDNA or cDNA will determine how the length of a feature is "
             "calculated for normalization. If cDNA, effective length will "
             "serve as feature length")
    parser.add_argument('--filter',
        dest='cat_only',
        action='store_true',
        help="only output abundances for features with an associated feature "
             "category [default: output all]")
    parser.add_argument('--nonunique',
        action='store_true',
        help="allow reads to align with more than one feature")
    compression = parser.add_mutually_exclusive_group()
    compression.add_argument('--gzip',
        dest='gzipped',
        action='store_true',
        help="compress output using the gzip algorithm")
    compression.add_argument('--bzip2',
        dest='bzipped',
        action='store_true',
        help="compress output using the bzip2 algorithm")
    compression.add_argument('--lzma',
        dest='lzma',
        action='store_true',
        help="compress output using the lzma algorithm")
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    # Argument sanity checks
    if (args.category and not args.map_files) or \
        (args.map_files and not args.category):
        parser.error("error: -m/--mapping and -c/--category must be supplied "
                     "together")

    if args.alignment_file == '-' and args.feature_file == '-':
        parser.error("error: standard input (stdin) can only be redirected to "
            "a single positional argument")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('count_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    if args.gzipped:
        compression = '.gz'
    elif args.bzipped:
        compression = '.bz2'
    elif args.lzma:
        compression = '.xz'
    else:
        compression = ''

    allowed_units = ["counts", "tpm", "custom", "prop", "fpk", "fpkm"]

    out_handles = {}
    for unit in args.norm:
        if unit not in allowed_units:
            print("warning: unknown metric of abundance '{}' provided to "
                  "-u/--unit. Please see the help message for a list of the "
                  "allowed units".format(unit), file=sys.stderr)
            continue

        outfile = "{}.{}.csv{}".format(args.outpref, unit, compression)
        try:
            out_h = open_io(outfile, mode='wb')
        except AttributeError:
            print("error: unable to write to '{}'".format(outfile), \
                  file=sys.stderr)
            sys.exit(1)

        out_handles[unit] = out_h

    if not out_handles:
        print("error: no output files can be created. Please re-run with one "
              "or more of the accepted units of abundance", file=sys.stderr)
        sys.exit(1)

    overlap_mode = args.mode
    minaqual = args.minqual
    feature_type = args.ftype
    id_field = args.attr
    category_field = args.category
    are_transcripts = args.transcripts
    category_only = args.cat_only
    multi_aln = args.nonunique

    match_types = ('M', '=', 'X')

    if args.aformat == "sam":
        align_reader = HTSeq.SAM_Reader
    else:  #must be BAM then
        align_reader = HTSeq.BAM_Reader

    if args.map_files:
        mapping = {}
        for map_file in args.map_files:
            json_map = json.load(open_io(map_file))
            mapping = {**json_map, **mapping}
    else:
        mapping = None

    # Store features in genomic arrays
    features = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    counts = {}

    # Iterate over GFF3 file, storing features to estimate coverage for
    no_attr = 0
    f_totals = 0
    ftype_totals = 0

    try:
        if args.feature_file == '-':
            gff = HTSeq.GFF_Reader(sys.stdin)
        else:
            gff = HTSeq.GFF_Reader(args.feature_file)

        for f in gff:
            f_totals += 1

            try:
                feature_id = f.attr[id_field]
            except KeyError:
                no_attr += 1
                feature_id = "unkwn_{:08}".format(no_attr)

            # Skip features of wrong type
            if feature_type:
                if f.type == feature_type:
                    ftype_totals += 1
                else:
                    continue

            # Store feature length for normalization
            feature_length = abs(f.iv.end - f.iv.start + 1)

            features[f.iv] += feature_id  #for mapping alignments
            counts[feature_id] = {'count': 0, 'length': feature_length}
    except:
        print("error: problem occured when processing GFF3 file at line {}"
              .format(gff.get_line_number_string()), file=sys.stderr)
        sys.exit(1)

    # Verify GFF3 file contains features of the specified type
    if ftype_totals == 0:
        print("error: no features of type '{}' found.\n"
            .format(feature_type), file=sys.stderr)
        sys.exit(1)

    if no_attr > 0:
        print("warning: found {!s} features without a '{}' attribute.\n"\
              .format(no_attr, id_field), file=sys.stderr)

    # Check alignment file formatting
    try:
        if args.alignment_file == '-':
            read_seq_file = align_reader(sys.stdin)
            read_seq_iter = iter(read_seq_file)
            first_read = next(read_seq_iter)
            read_seq = itertools.chain([first_read], read_seq_iter)
        else:
            read_seq_file = align_reader(args.alignment_file)
            read_seq = read_seq_file
            first_read = next(iter(read_seq))
    except:
        print("error: unable to read the alignment file. Please verify that "
              "the formatting is correct.", file=sys.stderr)
        sys.exit(1)

    pe_mode = first_read.paired_end  #reads are paired-end or single-end
    if pe_mode:
        if args.order == "name":
            read_seq = HTSeq.pair_SAM_alignments(read_seq)
        else:  #order is by position
            read_seq = HTSeq.pair_SAM_alignments_with_buffer(read_seq, \
                       max_buffer_size=args.buffer_size)

    # Iterate over alignment file
    empty = 0  #reads aligned somewhere in the assembly, but not to a feature
    duplicate = 0  #reads are duplicates of other reads
    ambiguous = 0  #reads overlapping more than one feature
    notaligned = 0  #unaligned reads
    lowqual = 0  #reads not passing minimum threshold for alignment quality
    nonunique = 0  #reads having multiple alignments with similar score
    r_totals = 0  #total reads
    aln_totals = 0  #correctly mapped to a feature
    fld = []  #fragment length / insert-size distribution 
    for r in read_seq:

        r_totals += 1

        if not pe_mode:  #single-end read mapping

            # Check if read aligned
            if not r.aligned:
                notaligned += 1
                continue

            # Check if the read aligned uniquely
            try:
                if r.optional_field("NH") > 1:
                    nonunique += 1
                    print("warning: read '{}' has multiple alignments with "
                          "similar score.\n".format(r.iv.chrom), \
                          file=sys.stderr)
                    continue
            except KeyError:
                pass

            # Cehck if the alignment passed the quality requirement
            if r.aQual < minaqual:
                lowqual += 1
                continue

            # Check whether the read was marked as a duplciate
            if r.pcr_or_optical_duplicate:
                duplicate += 1
                continue

            # Store read coordiantes
            iv_seq = (co.ref_iv for co in r.cigar if co.type in match_types \
                      and co.size > 0)

        else:  #paired-end read mapping

            # Store pair coordinates
            try:
                first_r, second_r = r
            except ValueError:
                notaligned += 1
                continue

            if first_r is None or second_r is None:
                notaligned += 1
                continue

            if first_r is not None and first_r.aligned:
                iv_seq = (co.ref_iv for co in first_r.cigar if co.type in \
                          match_types and co.size > 0)
            else:
                iv_seq = tuple()

            if second_r is not None and second_r.aligned:
                iv_seq = itertools.chain(iv_seq, (co.ref_iv for co in \
                         second_r.cigar if co.type in match_types and \
                         co.size > 0))
            else:
                if (first_r is None) or not (first_r.aligned):
                    notaligned += 1
                    continue

            # Check whether either read aligned more than once
            try:
                if (first_r.optional_field("NH") > 1) or \
                   (second_r.optional_field("NH") > 1):
                    nonunique += 1
                    print("warning: read '{}' has multiple alignments with "
                          "similar score.\n".format(first_r.iv.chrom), \
                          file=sys.stderr)
                    continue
            except KeyError:
                pass

            # Check if both reads passed the quality requirement
            if first_r.aQual < minaqual or second_r.aQual < minaqual:
                lowqual += 1
                continue

            # Check if the read pair was marked as a duplicate
            if first_r.pcr_or_optical_duplicate or \
                second_r.pcr_or_optical_duplicate:
                duplicate += 1
                continue

            # Append fragment length/insert-size to distribution
            try:
                fld.append(first_r.inferred_insert_size)
            except AttributeError:
                pass

        # Handle case where reads might overlap more than one feature
        try:
            if overlap_mode == "union":
                fs = set()  #store feature names when reads align
                for iv in iv_seq:
                     if iv.chrom not in features.chrom_vectors:
                         raise UnknownChrom

                     for iv2, fs2 in features[iv].steps():
                         fs = fs.union(fs2)

            else:  #intersection
                fs = None
                for iv in iv_seq:
                    if iv.chrom not in features.chrom_vectors:
                        raise UnknownChrom

                    for iv2, fs2 in features[iv].steps():
                        if len(fs2) > 0 or overlap_mode == "intersection-strict":
                            if fs is None:
                                fs = fs2.copy()
                            else:
                                fs = fs.intersection(fs2)
                                
            # If a read correctly mapped to a feature, increment its abundance
            if not fs:
                empty += 1
                continue
            elif len(fs) > 1:
                ambiguous += 1
                if not multi_aln:
                    continue
            else:
                aln_totals += 1

            for fsi in list(fs):
                counts[fsi]['count'] += 1

        except UnknownChrom:
            empty += 1

    unaln_totals = empty + ambiguous + lowqual + notaligned + nonunique + \
                      duplicate
    nmapped = aln_totals + empty + ambiguous + nonunique + lowqual

    for unit in out_handles:
        # Set scaling function
        if unit == 'fpk':
            norm_method = scale_abundance_fpk
            scaling_factor = None
        elif unit == 'fpkm':
            norm_method = scale_abundance_fpkm
            # Scaling factor is all mapped reads
            scaling_factor = nmapped
            print("info: the total number of mapped reads used in calculation "
                  "of FPKM is {!s}.\n".format(nmapped), file=sys.stderr)
        elif unit == 'tpm':
            norm_method = scale_abundance_tpm
            rates = [counts[j]['count'] / counts[j]['length'] for j in counts]
            rate_sum = sum(rates)
            print("info: the sum of all counts per bp rates used in "
                    "estimating fragment proportions is {:.2f}.\n"\
                  .format(rate_sum), file=sys.stderr)
            # Scaling factor is sum of all reads per base rates
            scaling_factor = rate_sum
        elif unit == 'custom':
            norm_method = scale_abundance_prop
            rates = [counts[j]['count'] / counts[j]['length'] for j in counts]
            rate_sum = sum(rates)
            print("info: the sum of all counts per bp rates used in "
                    "estimating fragment proportions is {:.2f}.\n"\
                  .format(rate_sum), file=sys.stderr)
            scaling_factor = args.sfactor / rate_sum
        elif unit == 'prop':
            norm_method = scale_abundance_prop
            rates = [counts[j]['count'] / counts[j]['length'] for j in counts]
            rate_sum = sum(rates)
            print("info: the sum of all counts per bp rates used in "
                    "estimating fragment proportions is {:.2f}.\n"\
                  .format(rate_sum), file=sys.stderr)
            scaling_factor = 1 / rate_sum
        else:  #default is counts
            norm_method = scale_abundance_none
            scaling_factor = None

        if are_transcripts and not pe_mode:
            print("warning: unable to calculate effective length from single-end "
                  "reads. Will use sequence length instead.\n", file=sys.stderr)
            calc_length = return_first_arg
        elif are_transcripts and pe_mode:
            calc_length = calculate_effective_length
        else:
            calc_length = return_first_arg

        out_h = out_handles[unit] 

        # Abundance normalization
        abundances = {}
        unkwn_feat = 0
        no_map = 0
        for feature in counts:

            fcount = counts[feature]['count']
            flen = calc_length(counts[feature]['length'], fld)

            feature_abundance = norm_method(fcount, flen, scaling_factor)

            # Map to higher order features, if applicable
            if category_field:
                try:
                    # Ensure that feature has corresponding entry in database
                    feature_map = mapping[feature]
                except KeyError:
                    no_map += 1
                    if not category_only:
                        # Keep all features, even the uncategorized ones
                        abundances[feature] = abundances.get(feature, 0) + \
                                              feature_abundance
                    continue
                else:
                    try:
                        # Ensure that entry has relevant category field
                        category = feature_map[category_field]
                    except KeyError:
                        unkwn_feat += 1
                        if not category_only:
                            abundances[feature] = abundances.get(feature, 0) + \
                                                  feature_abundance
                        continue

                # Handle case where feature has more than one category, such 
                # as if a protein sequence is assigned to more than one gene 
                # family
                for category in category.split(','):
                    abundances[category.lstrip()] = \
                        abundances.get(category, 0) + feature_abundance

            else:
                abundances[feature] = abundances.get(feature, 0) + \
                                      feature_abundance

        # "UNMAPPED" can be interpreted as a single unknown gene of length one
        # kilobase recruiting all reads that failed to map to input features
        #abundances['UNMAPPED'] = unaln_totals

        # Output abundances sorted by key name
        for fn in sorted(abundances):
            if not fn.startswith("unkwn_"):
                write_io(out_h, "{}\t{!s}\n".format(fn, abundances[fn]))

        out_h.close()

    if unkwn_feat > 0:
        print("warning: found '{!s}' features without the '{}' field in the "
              "relational database.\n".format(unkwn_feat, category_field), \
              file=sys.stderr)

    if no_map > 0:
        print("warning: found {!s} features without an entry in the "
              "relational database.\n".format(no_map), file=sys.stderr)

    # Output statistics
    print("Features processed:", file=sys.stderr)
    print("  - feature totals:\t{!s}".format(f_totals), file=sys.stderr)
    if feature_type:
        print("  - of relevant type:\t{!s}".format(ftype_totals), \
              file=sys.stderr)
    print("  - unique features:\t{!s}".format(len(counts)), file=sys.stderr)
    print("Reads processed:", file=sys.stderr)
    print("  - read totals:\t{!s}".format(r_totals), file=sys.stderr)
    print("  - successfully mapped:\t{!s}".format(aln_totals), \
          file=sys.stderr)
    if multi_aln:
        print("    - ambiguous alignment:\t{!s}".format(ambiguous), \
              file=sys.stderr)
    print("  - unsuccessfully mapped:\t{!s}".format(unaln_totals), \
          file=sys.stderr)
    print("    - no feature\t{!s}".format(empty), file=sys.stderr)
    if not multi_aln:
        print("    - ambiguous alignment\t{!s}".format(ambiguous), \
              file=sys.stderr)
    print("    - too low alignment quality\t{!s}".format(lowqual), \
          file=sys.stderr)
    print("    - not aligned\t{!s}".format(notaligned), file=sys.stderr)
    print("    - duplicate\t{!s}".format(duplicate), file=sys.stderr)
    print("    - alignment not unique\t{!s}".format(nonunique), \
          file=sys.stderr)
    print("", file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to count {!s} fragments for {!s} features"\
          .format(total_time, r_totals, f_totals), file=sys.stderr)
    print("", file=sys.stderr)


if __name__ == "__main__":
   main()
   sys.exit(0)
