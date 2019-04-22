#! /usr/bin/env python
"""Functions for handling input files
"""

from bio_utils.iterators import B6Reader, B6Entry
from bz2 import BZ2File
from gzip import GzipFile
from lzma import LZMAFile
import io
import os
import sys


class FormatError(Exception):
    """A simple exception that is raised when an input object is formatted
    incorrectly
    """
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class InputError(Exception):
    """A simple exception that is raised when bad argument values encountered
    """
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


class Homolog:
    """A class to store homology search results. Provides support for the 
    output formats of different search methods, including the B6 format of 
    BLAST+ and the tblout format of the protein search programs from HMMER

    Attributes:
        seqid (str): identifier of the dataset sequence. Typically a homology 
            search method will query a sequence or collection of sequences 
            against a database of sequences, profile HMMs, PSSMs, etc; 
            however, some methods (such as hmmsearch) reverse the role of 
            query and subject. For consistency this attribute should contain 
            the identifier of the dataset sequence, regardless of role assumed 
            in the search

        dbid (str): identifier of the database entry

        evalue (float): expectation value for the target/query comparsion

        score (float): bit score for the target/query comparsion

        seq_start (int): alignment start position for dataset sequence

        seq_end (int): alignment end position for dataset sequence

        db_start (int): alignment start position for database entry

        db_end (int): alignment end position for database entry
    """

    def __init__(self):
        """Initialize variables to store entry data"""
        self.seqid = None
        self.dbid = None
        self.evalue = None
        self.score = None
        self.seq_start = None
        self.seq_end = None
        self.db_start = None
        self.db_end = None
        self.description = None

    def parse_entry(self, entry, reverse_query: bool=False, domain: bool=False):
        """Convert B6Entry or HMMEntry objects to Homolog object. 

        Args:
            entry (class): object containing the results of a single 
                query/target comparison. Can be in the form of either a 
                B6Entry or HHMEntry object

            revere_query (bool): dataset sequence assumed the role of target 
                in homology search if True. Only applies to HMMEntry objects

            domain (bool): use domain results if True, else use sequence 
                results. Only applies to HMMEntry objects
        """
        if isinstance(entry, B6Entry):
            self.seqid = entry.query
            self.dbid = entry.subject
            self.evalue = entry.evalue
            self.score = entry.bitscore
            self.seq_start = entry.query_start
            self.seq_end = entry.query_end
            self.db_start = entry.subject_start
            self.db_end = entry.subject_end

        elif isinstance(entry, HMMEntry):
            if reverse_query:
                self.seqid = entry.tacc if entry.tacc else entry.target
                self.dbid = entry.qacc if entry.qacc else entry.query
            else:
                self.seqid = entry.qacc if entry.qacc else entry.query
                self.dbid = entry.tacc if entry.tacc else entry.target
            if domain:
                self.evalue = entry.dom_evalue
                self.score = entry.dom_score
            else:
                self.evalue = entry.seq_evalue
                self.score = entry.seq_score

        else:
            raise FormatError("unsupported object class")


class HMMEntry:
    """A class to store alignments from the target hits table of HMMER, which 
    consists of a single line for each query / target comparison passing the
    reporting threshold

    Attributes:
        target (str): target name

        tacc (str): target accession

        query (str): query name

        qacc (str): query accession

        seq_evalue (float): e-value of full sequence

        seq_score (int): corrected bit score for target / query comparison

        seq_bias (int): bias correction applied to bit score of full sequence

        dom_evalue (float): e-value of best-scoring domain

        dom_score (int): corrected bit score of the best-scoring domain

        dom_bias (int): bias correction applied to bit score of best-scoring 
            domain

        reg (tuple): tuple with two items corresponding to fields reg and clu 
            in the hits table. These fields are defined as number of discrete 
            regions identified and the number of regions that appear to be 
            mulitdomain

        env (tuple): tuple with two items corresponding to the fields ov and 
            env in the hits table. These fields are defined as the number of 
            envelopes defined by stochastic traceback clustering that overlap 
            other envelopes and the total number of envelopes identified

        dom (tuple): tuple with four items corresponding to fields exp, dom, 
            rep and inc in the hits table. These fields are defined as the 
            expected number of domains, total number of domains defined, 
            number of domains passing the reporting threshold, and number 
            of domains satisfying the inclusion threshold

        description (str): description of the target sequence
    """

    def __init__(self):
        """Initialize variables to store entry data"""
        self.target = None
        self.tacc = None
        self.query = None
        self.qacc = None
        self.seq_evalue = None
        self.seq_score = None
        self.seq_bias = None
        self.dom_evalue = None
        self.dom_score = None
        self.dome_bias = None
        self.reg = None
        self.env = None
        self.dom = None
        self.description = None

    def write(self, *args):
        """Restore entry to original format (with the exception that output 
        is tab-delimited instead of justified)

        Returns:
            str: properly formatted string containing the target / query 
                comparison
        """
        # Format entry for writing
        reg, clu = self.reg
        ov, env = self.env
        exp, dom, rep, inc = self.dom
        fields = [self.target, self.tacc, self.query, self.qacc, \
            self.seq_evalue, self.seq_score, self.seq_bias, self.dom_evalue, \
            self.dom_score, self.dom_bias, exp, reg, clu, ov, env, dom, rep, \
            inc, self.description]

        fstr = "\t".join(fields)
        return '{}{}'.format(fstr, os.linesep)


class HMMReader():
    """Class to read from HMMER tblout files and store lines as HMMEntry 
    objects

    Attributes:
        handle (file): tblout file handle, can be any iterator so long as it
            it returns subsequent target / query comparisons (in the form of
            a new line)

        filename (str): name of the file
    
        current_line (int): current line in file [default: 0]
    """

    def __init__(self, handle):
        """Initialize variables to store file information"""

        self.handle = handle
        self.filename = handle.name
        self.current_line = 0

    def iterate(self, start_line=None, comments: bool=False):
        """Iterate over HMMER tblout file and return single hits

        Args:
            start_line (str): Next HMM entry. If 'handle' has been partially
                read and you want to start iterating at the next entry, read 
                the next HMM entry and pass it to this variable when calling 
                the method

            comments (bool): Yields comments if True, else skips lines starting
                with "#"

        Yields:
            HMMEntry: class containing a single target / query comparison
        """

        handle = self.handle

        # Speed tricks: reduces function calls
        split = str.split
        strip = str.strip

        # Begin reading text
        if start_line is None:
            line = next(handle)  # Read first HMM entry
        else:
            line = start_line  # Set header to given line

        # Check if input is text or bytestream
        if (isinstance(line, bytes)):
            def next_line(i):
                return next(i).decode('utf-8')

            line = strip(line.decode('utf-8'))
        else:
            next_line = next
            line = strip(line)

        try:  # Manually construct a for loop to improve speed by using 'next'

            while True:  # Loop until StopIteration Exception raised

                self.current_line += 1

                data = HMMEntry()

                if line.startswith('#') and not comments:
                    line = strip(next_line(handle))
                    continue
                elif line.startswith('#') and comments:
                    yield line
                    line = strip(next_line(handle))
                    continue

                split_line = split(' '.join(line.split()), ' ', maxsplit=18)

                # Replace empty values with None
                fields = [None if i == '-' else i for i in split_line]

                data.target = fields[0]
                data.tacc = fields[1]
                data.query = fields[2]
                data.qacc = fields[3]
                data.seq_evalue = fields[4]
                data.seq_score = fields[5]
                data.seq_bias = fields[6]
                data.dom_evalue = fields[7]
                data.dom_score = fields[8]
                data.dom_bias = fields[9]
                data.reg = (fields[11], fields[12])
                data.env = (fields[13], fields[14])
                data.dom = (fields[10], fields[15], fields[16], fields[17])
                data.description = fields[18]

                line = strip(next_line(handle))  # Raises StopIteration at EOF

                yield data

        except StopIteration:  # Yield last entry
            if data.target:
                yield data
            else:  #handle case where file ends in comment
                pass


def open_io(infile, **kwargs):
    """Open input files based on file extension for reading or writing
    """

    algo = io.open  # Default to plaintext

    algo_map = {
                'bz2': BZ2File,
                'gz': GzipFile,
                'xz': LZMAFile
               }

    # Base compression algorithm on file extension
    ext = infile.split('.')[-1]
    try:
        algo = algo_map[ext]
    except KeyError:
        pass

    handle = algo(infile, **kwargs)

    return handle

def write_io(out_h, output):
    """Write output to file handle
    """
    try:  # Output to file
        out_h.write(output.encode('utf-8'))
    except TypeError:  # Output to stdout
        out_h.write(output)

