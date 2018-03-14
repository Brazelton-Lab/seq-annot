#!/usr/bin/env python
"""
Calculate genomic feature abundances.

Usage:
    count_features [options] [-o out.csv] [-m in.json] in.aln in.gff3

Required inputs are a GFF3 file of annotated features and an alignments file 
in SAM or BAM format. Optional input is a relational database containing 
features mapped to feature categories, such as genes to gene families. Output 
is a tabular file of feature abundances.

The compression algorithm is automatically detected for input files based on 
the file extension. To compress output, add the appropriate file extension 
to the output file name (e.g. .gz, .bz2). Leave off '--out' to direct output 
to standard output (stdout).

Copyright:

    count_features Counts coverage of GFF file from SAM/BAM file
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
from seq_annot.seqio import open_input
from statistics import mean
import sys
import textwrap
from time import time

__author__ = 'Christopher Thornton'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Beta"
__version__ = '1.3.0'


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


def scale_abundance_biomass(counts, length, sfactor, *args):
    """Normalize coverage by biomass and feature length.

    formula:
        biomass_i = countPerBase_i * ( numCells/totalCountsPerBase )

    countsPerBase - number of counts per base ( counts/effLen )
    counts - number of reads mapped to a feature (i.e. a gene sequence)
    effLen - effective length of a feature
    totalCountsPerBase - sum of all counts per base rates
    numCells - estimated number of cells in the sample
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
        help="input feature annotation file in GFF3 format")
    parser.add_argument('-m', '--mapping', 
        metavar='in.json',
        dest='map_files',
        action=ParseSeparator,
        sep=',',
        help="input relational database in JSON format containing features "
             "mapped to feature categories, such as genes to gene families "
             "or exons to genes. Abundance estimates for feature categories "
             "will be reported instead of for features")
    parser.add_argument('-o', '--out',
        metavar='out.csv',
        action=Open,
        mode='wt',
        default=sys.stdout,
        help="output tabular file of feature abundances [default: output "
             "to stdout]")
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
        help="alignment file sorting scheme [default: position]. The input "
             "alignments file must be pre-sorted either by position or by "
             "read name. This will be ignored for single-end reads. Options "
             "are 'position' or 'name'" )
    parser.add_argument('-t', '--type', 
        metavar='TYPE', 
        dest='ftype',
        default='CDS',
        help="feature type (3rd column in GFF file) to estimate abundance for "
             "[default: CDS]. All features of other type will be ignored")
    parser.add_argument('-a', '--attr',
        metavar='ATTRIBUTE',
        default="Name",
        help="GFF attribute to use as the ID for the calculated abundance "
             "[default: 'Name']. The value will also be used as the search "
             "ID for the relational database")
    parser.add_argument('-e', '--mode', 
        metavar='MODE',
        choices=["union", "intersection-strict", "intersection-nonempty"],
        default="union",
        help="mode to handle reads that overlap more than one feature "
            "[default: union]. Options are 'union', 'intersection-strict', "
            "and 'intersection-nonempty'.")
    parser.add_argument('-n', '--norm', 
        metavar='METHOD',
        choices=['none', 'fpk', 'fpkm', 'tpm', 'biomass'], 
        default='none',
        help="normalization method to use in abundance estimation [default: "
             "none]. Choices are 'none', 'fpk' (fragments per kilobase of "
             "feature), 'fpkm' (fragements per kilobase of feature per "
             "million mapped fragments), 'tpm' (transcripts per million), and "
             "'biomass'. For fpk, fpkm, and tpm, feature length will be taken "
             "from the relational database if provided, otherwise it will be "
             "calculated from feature start and end fields in the GFF3 file. "
             "FPK is a normalization method to correct for differences in "
             "feature lengths within a sample. It should only be used to "
             "compare features within a single sample. FPKM and TPM are "
             "normalization methods for correcting differences in both "
             "feature lengths within a sample and for sequencing effort "
             "(sample size) between samples. The functional difference "
             "between TPM and FPKM is that TMP is a proportional measurement, "
             "making it easier to identify the extent that the relative "
             "'importance' of a given feature changes between samples")
    parser.add_argument('-g', '--nanograms',
        metavar='TOTAL',
        dest='total_dna',
        type=float,
        help="sample total DNA in nanograms. Required when the specified "
             "normalization method is 'biomass'")
    parser.add_argument('-d', '--dna',
        metavar='AVG',
        dest='avg_dna',
        type=float,
        default=0.000002,
        help="average concentration of DNA per cell [default: 2E-6 (ng/cell)]. "
             "Should be used in conjunction with --norm biomass")
    parser.add_argument('-q', '--qual', 
        metavar='THRESHOLD',
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
    parser.add_argument('--version',
        action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if args.category and not args.map_files:
        parser.error("error: -m/--mapping must also be supplied when "
                     "-c/--category is used")

    if args.norm == "biomass" and not args.total_dna:
        parser.error("error: -g/--nanograms is required when "
                     "normalizing by biomass")

    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('count_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)

    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs
    out_h = args.out.write

    overlap_mode = args.mode
    minaqual = args.minqual
    feature_type = args.ftype
    id_field = args.attr
    category_field = args.category
    are_transcripts = args.transcripts
    category_only = args.cat_only

    if args.aformat == "sam":
        align_reader = HTSeq.SAM_Reader
    else:  #must be BAM then
        align_reader = HTSeq.BAM_Reader

    if args.map_files:
        mapping = {}
        for map_file in args.map_files:
            json_map = json.load(open_input(map_file))
            mapping = {**json_map, **mapping}
    else:
        mapping = None

    # Iterate over GFF3 file, storing feature counts into a dictionary
    features = HTSeq.GenomicArrayOfSets("auto", False)
    counts = {}

    no_attr = 0
    feature_totals = 0
    feature_type_totals = 0
    try:
        gff = HTSeq.GFF_Reader(args.feature_file)

        for f in gff:
            feature_totals += 1

            try:
                feature_id = f.attr[id_field]
            except KeyError:
                no_attr += 1
                feature_id = "unkwn_{:08}".format(no_attr)

            # Skip features of other type
            if feature_type:
                if f.type == feature_type:
                    feature_type_totals += 1
                else:
                    continue

            # Store feature length for normalization
            feature_length = abs(f.iv.end - f.iv.start)

            features[f.iv] += feature_id  #for mapping alignments
            counts[feature_id] = {'count': 0, 'length': feature_length}
    except:
        print("error: problem occured when processing GFF3 file at line ({})"
              .format(gff.get_line_number_string()), file=sys.stderr)
        sys.exit(1)

    # Verify GFF3 file contains features of the specified type
    if feature_type_totals == 0:
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
    ambiguous = 0  #reads overlapping more than one feature
    notaligned = 0  #unaligned reads
    lowqual = 0  #reads not passing minimum threshold for alignment quality
    nonunique = 0  #reads having multiple alignments with similar score
    aln_totals = 0  #total reads
    mapped_totals = 0  #correctly mapped to a feature
    fld = []  #fragment length / insert-size distribution 
    for r in read_seq:

        aln_totals += 1

        if not pe_mode:  #single-end read mapping

            # Did the read align at all?
            if not r.aligned:
                notaligned += 1
                continue

            # Did the read align uniquely?
            try:
                if r.optional_field("NH") > 1:
                    nonunique += 1
                    print("warning: read '{}' has multiple alignments with "
                          "similar score. It is recommended that best-hit "
                          "filtering is performed prior to abundance "
                          "estimation.\n".format(r.iv.chrom), file=sys.stderr)
                    continue
            except KeyError:
                pass

            # Did the alignment pass the minimum quality threshold?
            if r.aQual < minaqual:
                lowqual += 1
                continue

            iv_seq = (invert_strand(co.ref_iv) for co in r.cigar if \
                      co.type == "M" and co.size > 0)

        else:  #paired-end read mapping
            # Did both reads properly align?
            if r[0] is not None and r[0].aligned:
                iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar if \
                          co.type == "M" and co.size > 0)
            else:
                iv_seq = tuple()

            if r[1] is not None and r[1].aligned:
                iv_seq = itertools.chain(iv_seq, (co.ref_iv for co in \
                         r[1].cigar if co.type == "M" and co.size > 0))
            else:
                if (r[0] is None) or not (r[0].aligned):
                    notaligned += 1
                    continue

            # Did either one of the reads align more than once?
            try:
                if (r[0] is not None and r[0].optional_field("NH") > 1 ) or \
                   (r[1] is not None and r[1].optional_field("NH") > 1):
                    nonunique += 1
                    print("warning: read '{}' has multiple alignments with "
                          "similar score. It is recommended that best-hit "
                          "filtering is performed prior to abundance "
                          "estimation.\n".format(r[0].iv.chrom), \
                          file=sys.stderr)
                    continue
            except KeyError:
                pass

            # Did either one of the read pair align poorly?
            if (r[0] and r[0].aQual < minaqual) or (r[1] and r[1].aQual < \
                minaqual):
                lowqual += 1
                continue

        # Handle case where reads might overlap more than one feature
        try:
            if overlap_mode == "union":
                fs = set()
                for iv in iv_seq:
                     if iv.chrom not in features.chrom_vectors:
                         raise UnknownChrom

                     for iv2, fs2 in features[iv].steps():
                         fs = fs.union(fs2)

            elif overlap_mode == "intersection-strict" or \
                overlap_mode == "intersection-nonempty":
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
            elif len(fs) > 1:
                ambiguous += 1
            else:
                counts[list(fs)[0]]['count'] += 1
                mapped_totals += 1

                #append fragment length/insert-size to fld
                try:  #paired-end reads
                    fld.append(r[0].inferred_insert_size)
                except AttributeError:  #no mate
                    continue
                    #fld.append(len(r[0].read.seq))

        except UnknownChrom:
            empty += 1

    unmapped_totals = empty + ambiguous + lowqual + notaligned + nonunique
    nmapped = mapped_totals + empty + ambiguous + nonunique + lowqual

    # Set scaling function
    if args.norm == 'fpk':
        norm_method = scale_abundance_fpk
        scaling_factor = None
    elif args.norm == 'fpkm':
        norm_method = scale_abundance_fpkm
        # Scaling factor is all mapped reads
        scaling_factor = nmapped
    elif args.norm == 'tpm':
        norm_method = scale_abundance_tpm
        rates = [counts[i]['count'] / counts[i]['length'] for i in counts]
        print(sum(rates), file=sys.stderr)
        # Scaling factor is sum of all reads per base rates
        scaling_factor = sum(rates)
    elif args.norm == 'biomass':
        norm_method = scale_abundance_biomass
        num_cells = args.total_dna / args.avg_dna
        rates = [counts[i]['count'] / counts[i]['length'] for i in counts]
        # Scaling factor is number cells over sum of reads per base rates
        scaling_factor = num_cells / sum(rates)
    else:  #default is none
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

    # Abundance normalization
    abundances = {}
    unkwn_feat = 0
    no_map = 0
    for feature in counts:

        fcount = counts[feature]['count']
        efflen = calc_length(counts[feature]['length'], fld)

        feature_abundance = norm_method(fcount, efflen, scaling_factor)

        # Map to higher order features, if applicable
        if category_field:
            try:
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
                    category = feature_map[category_field]
                except KeyError:
                    unkwn_feat += 1
                    if not category_only:
                        abundances[feature] = abundances.get(feature, 0) + \
                                              feature_abundance
                    continue

            # Handle if feature has more than one category
            for category in category.split(','):
                abundances[category.lstrip()] = abundances.get(category, 0) + \
                                                feature_abundance

        else:
            abundances[feature] = abundances.get(feature, 0) + feature_abundance

    if unkwn_feat > 0:
        print("warning: found '{!s}' features without the '{}' field in the "
              "relational database.\n".format(unkwn_feat, category_field), \
              file=sys.stderr)

    if no_map > 0:
        print("warning: found {!s} features without an entry in the "
              "relational database.\n".format(no_map), file=sys.stderr)

    # "UNMAPPED" can be interpreted as a single unknown gene of length one
    # kilobase recruiting all reads that failed to map to input features
    #abundances['UNMAPPED'] = unmapped_totals

    # Output abundances
    for fn in sorted(abundances):
        if not fn.startswith("unkwn_"):
            out_h("{}\t{!s}\n".format(fn, abundances[fn]))

    # Output statistics
    print("Total number of features:\t{!s}".format(feature_totals),\
          file=sys.stderr)
    if feature_type:
        print("  - features of relevant type:\t{!s}"\
              .format(feature_type_totals), file=sys.stderr)
    print("Total number of reads:\t{!s}".format(aln_totals),\
          file=sys.stderr)
    print("  - number of mapped reads:\t{!s}".format(nmapped),\
          file=sys.stderr)
    print("  - number of reads that successfully mapped to a feature:\t{!s}"\
          .format(mapped_totals), file=sys.stderr)
    print("  - number of reads that failed to correctly map to any feature:\t{!s}"\
          .format(unmapped_totals), file=sys.stderr)
    print("    - no feature\t{!s}".format(empty), file=sys.stderr)
    print("    - ambiguous alignment\t{!s}".format(ambiguous), file=sys.stderr)
    print("    - too low alignment quality\t{!s}".format(lowqual), \
          file=sys.stderr)
    print("    - not aligned\t{!s}".format(notaligned), file=sys.stderr)
    print("    - alignment not unique\t{!s}\n".format(nonunique), \
          file=sys.stderr)

    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to count {!s} fragments for {!s} features\n"\
          .format(total_time, aln_totals, feature_totals), file=sys.stderr)


if __name__ == "__main__":
   main()
   sys.exit(0)
