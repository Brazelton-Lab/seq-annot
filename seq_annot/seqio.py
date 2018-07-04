#! /usr/bin/env python
"""Functions for handling input files
"""

from bz2 import BZ2File
from gzip import GzipFile
from lzma import LZMAFile
import io
import sys

def open_io(infile, mode='rb'):
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

    handle = algo(infile, mode=mode)

    return handle
