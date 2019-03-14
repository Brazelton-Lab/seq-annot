#! /usr/bin/env python
"""Functions for handling input files
"""

from bz2 import BZ2File
from gzip import GzipFile
from lzma import LZMAFile
import io
import sys


class FormatError(Exception):
    """A simple exception that is raised when an input file is formatted
    incorrectly
    """
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)


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

