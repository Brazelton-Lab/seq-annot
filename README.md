# seq-annot
Tools to facilitate annotation and comparison of genomes and metagenomes

## About

seq-annot is a python package for annotating and counting genomic features 
in genomes and metagenomes. 

# Licence:

seq-annot is under [open source licence GPLv3](https://opensource.org/licenses/GPL-3.0)

## Requirements

Python 3.4+

Python Libraries:

* arandomness
* bio_utils
* HTSeq
* numpy

## Installation

pip install seq-annot

## Usage

### reldb
Integral to the annotation process is a sequence database with high-quality 
supporting information. seq-annot will accept supporting data in the form of a 
JSON-formatted relational database that provides context for entries in the 
sequence database, such as function, alternative IDs, organism, etc. Database 
entries can contain any number of fields. reldb is a tool used to perform 
various manipulations on one or more relational databases, including: entry 
merging, database subsetting, and format conversion.

#### Examples

    Find all entries with istA in the product description and output results as \
    a tabular CSV:
        reldb --csv --case -s "product:istA" in.mapping

    Merge entries by value of the field 'gene':
        reldb -m gene in.mapping

### screen_features
screen_features is for screening the results of a homology search using 
alternative phenotype-conferring snps, bitscore thresholds, and other scoring 
metrics.

#### Examples

    Filter hits by bitscore threshold specified in the relational database:
        screen_features -m in.mapping -b bitscore in.b6

### annotate_features
annotate_features is for annotating a GFF3 file of putative protein-coding 
genes using the results of a homology search to a database of genes with known 
or predicted function.

#### Examples

    Annotate predicted protein-coding genes, adding fields Dbxref and gene 
    from the relational database as additional attributes:
        annotate_features -m mapping.json -f Dbxref,gene -b in.b6 in.gff

### combine_features
combine_features is for combining GFF3 files containing feature annotations. 
Overlap conflicts can be resolved through a hierarchy of precedence determined 
by file input order.

#### Examples

    Merge two GFF3 files, using input file order to resolving overlap conflicts
    but allowing overlaps if found on different strands:
        combine_features --precedence --stranded in.gff1 in.gff2

### count_features
count_features is for estimating the abundance of genomic features using a 
GFF3 file of feature annotations and an read alignment file in SAM/BAM format.

#### Examples

    Estimate abundance for features in GFF in units of FPKM:
        count_features --nonunique -u fpkm in.bam in.gff

### compare_features
compare_features is for comparing annotated features across multiple samples.

#### Examples

    Create table of feature abundances by sample:
        compare_features -n sample1,sample2 in.csv1 in.csv2
    
### colocate_features
colocate_features is for finding evidence of co-residence between two sets of 
genomic features in an assembly of mixed populations (such as a metagenome 
assembly).

#### Examples

    Find all instances of co-residence between sets of features:
        colocate_features -a Alias -f in.fasta in.sets in.gff

