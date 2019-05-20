#! /usr/bin/env python
"""Custom ArgParse actions
"""

import argparse
from importlib import import_module
from inspect import getfullargspec
import io
from os import linesep
import sys
from warnings import warn

__author__ = 'Alex Hyer'
__email__ = 'theonehyer@gmail.com'
__license__ = 'GPLv3'
__maintainer__ = 'Christopher Thornton'
__status__ = 'Production/Stable'
__version__ = '1.2.0'


class ParseSeparator(argparse.Action):
    """Argparse Action that parses arguments by a user-supplied separator
       character

    Attributes:
        option_strings (list): list of str giving command line flags that
            call this action

        dest (str): Namespace reference to value

        sep (str): arbitrary separator character

        nargs (str): number of args as special char or int

        **kwargs (various): optional arguments to pass to super call
    """


    def __init__(self, option_strings, dest, sep=',', nargs=None, **kwargs):
        """Initialize class and spawn self as Base Class w/o nargs

        This class will "make" nargs by parsing the separator so it only
        accepts a single string, not a list.

        Args:
            option_strings (list): list of str giving command line flags that
                call this action

            dest (str): Namespace reference to value

            sep (str): arbitrary separator character

            nargs (str): number of args as special char or int

            **kwargs (various): optional arguments to pass to super call
        """

        # Only accept a single value to analyze
        if nargs is not None:
            raise ValueError('nargs not allowed for ParseSeparator')

        # Call self again but without nargs
        super(ParseSeparator, self).__init__(option_strings, dest, **kwargs)

        # Store and establish variables used in __call__
        self.kwargs = kwargs
        self.sep = sep

    def __call__(self, parser, namespace, value, option_string=None, **kwargs):
        """Called by Argparse when user specifies a character-separated list

        Simply splits a list by the supplied character and adds the values to
        namespace.

        Args:
            parser (ArgumentParser): parser used to generate values

            namespace (Namespace): namespace to set values for

            value (str): actual value specified by user

            option_string (str): argument flag used to call this function

        Raises:
            TypeError: if value is not a string

            ValueError: if value cannot, for any reason, be parsed
                by the separator character
        """

        # This try/except should already be taken care of by Argparse
        try:
            assert type(value) is str
        except AssertionError:
            raise TypeError('{0} is not a string'.format(value))

        try:
            arguments = [i.lstrip() for i in value.split(self.sep)]
        except:
            raise ValueError('{0} could not be parsed by {1}'
                             .format(value, self.sep))

        setattr(namespace, self.dest, arguments)


class Open(argparse.Action):
    """Argparse Action that detects and opens compressed files for rw

    Attributes:
        option_strings (list): list of str giving command line flags that
            call this action

        dest (str): Namespace reference to value

        mode (str): mode to pass to (de)compression algorithm

        nargs (bool): True if multiple arguments specified

        **kwargs (various): optional arguments to pass to argparse and algo
    """

    def __init__(self, option_strings, dest, mode='rb', nargs=None, **kwargs):
        """Initialize class and spawn self as Base Class w/o nargs

        Warns:
            ImportError: if Open cannot import a compression library, it warns 
                the user that it cannot open the corresponding file type

        Raises:
            ValueError: if nargs is not None, Open does not accept nargs
        """

        # Only accept a single value to analyze
        if nargs is not None:
            raise ValueError('nargs not allowed for Open')

        # Call self again but without nargs
        super(Open, self).__init__(option_strings, dest, **kwargs)

        # Store and establish variables used in __call__
        self.kwargs = kwargs
        self.mode = mode.lower().strip()
        self.modules = {}

        modules_to_import = {
            'bz2': 'BZ2File',
            'gzip': 'GzipFile',
            'lzma': 'LZMAFile'
        }

        # Dynamically import compression libraries and warn about failures
        for mod, _class in modules_to_import.items():
            try:
                self.modules[_class] = getattr(import_module(mod), _class)
            except (ImportError, AttributeError) as e:
                self.modules[_class] = open
                warn('Cannot process {0} files due to following error:'
                     '{1}{2}{1}You will need to install the {0} library to '
                     'properly use these files. Currently, such files will '
                     'open in text mode.'.format(mod, linesep, e))

    # Credits: https://stackoverflow.com/questions/13044562/
    # python-mechanism-to-identify-compressed-file-type-and-uncompress
    def __call__(self, parser, namespace, value, option_string=None, **kwargs):
        """Detects and opens compressed files

        Args:
            parser (ArgumentParser): parser used to generate values

            namespace (Namespace): namespace to set values for

            value (str): actual value specified by user

            option_string (str): argument flag used to call this function

            **kwargs (various): optional arguments later passed to the
                compression algorithm
        """

        filename = value  # For readability

        algo = io.open  # Default to plaintext

        # Capture any mode that isn't read, such as write or append
        if self.mode.lstrip('U')[0] != 'r':

            algo_map = {
                'bz2': self.modules['BZ2File'],
                'gz':  self.modules['GzipFile'],
                'xz':  self.modules['LZMAFile']
            }

            # Base compression algorithm on file extension
            ext = value.split('.')[-1]
            try:
                algo = algo_map[ext]
            except KeyError:
                pass

        # Basically read mode
        else:

            file_sigs = {
                b'\x42\x5a\x68': self.modules['BZ2File'],
                b'\x1f\x8b\x08': self.modules['GzipFile'],
                b'\xfd7zXZ\x00': self.modules['LZMAFile']
                }

            max_len = max(len(x) for x in file_sigs.keys())

            # Check if file is connected to a tty device
            if sys.stdin.isatty():
                # Check beginning of file for signature
                with io.open(filename, 'rb') as in_handle:
                    start = in_handle.read(max_len)
                    for sig in file_sigs.keys():
                        if start.startswith(sig):
                            algo = file_sigs[sig]
                            break


        # Filter all **kwargs by the args accepted by the compression algo
        algo_args = set(getfullargspec(algo).args)
        good_args = set(self.kwargs.keys()).intersection(algo_args)
        _kwargs = {arg: self.kwargs[arg] for arg in good_args}


        # Open the file using parameters defined above and store in namespace
        try:
            handle = algo(value, mode=self.mode, **_kwargs)
        except ValueError:
            mode = self.mode.lstrip('U')[0]
            handle = io.TextIOWrapper(algo(value, mode=mode, **_kwargs))
        except FileNotFoundError:
            warn('attempting to open {0} using {1} in mode {2}'\
                .format(value, algo, self.mode))

        setattr(namespace, self.dest, handle)
