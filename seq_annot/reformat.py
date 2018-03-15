#! /usr/bin/env python
"""

Copyright:

    reformat_features convert predicted features between formats
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
from bio_utils.iterators import b6_iter, gff3_iter
import json
from seq_annot.seqio import open_input
import sys
import textwrap
from time import time

__author__ = "Christopher Thornton"
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = "Alpha"
__version__ = "0.0.1"

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    args = parser.parse_args()


    # Output run information
    all_args = sys.argv[1:]
    print("{} {!s}".format('annotate_features', __version__), file=sys.stderr)
    print(textwrap.fill("Command line parameters: {}"\
          .format(' '.join(all_args)), 79), file=sys.stderr)
    print("", file=sys.stderr)


    # Track program run-time
    start_time = time()

    # Assign variables based on user inputs


    # Calculate and print statistics
    print("Total matches processed:\t{!s}".format(hits_totals), file=sys.stderr)
    print("Total features processed:\t{!s}".format(gff_totals), file=sys.stderr)
    print("  - features matching a single reference:\t{!s}".format(passed_totals), \
          file=sys.stderr)
    print("  - features with more than one matching reference:\t{!s}\n"\
          .format(conflict_totals), file=sys.stderr)


    # Calculate and print program run-time
    end_time = time()
    total_time = (end_time - start_time) / 60.0
    print("It took {:.2e} minutes to annotate {!s} features\n"\
          .format(total_time, gff_totals), file=sys.stderr)


if __name__ == "__main__":
    main()
    sys.exit(0)
