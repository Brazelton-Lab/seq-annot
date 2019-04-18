#! /usr/bin/env python
"""Functions for manipulating the relational databases
"""

import json
import re
from seq_annot.seqio import open_io, InputError, FormatError
import sys

def load_dbs(infiles:list, fields:list=None, csv:bool=False):
    """Load relational databases into memory

    Args:
        infiles (list): list of input relational databases in JSON format

        fields (list): list of fields to load [default: load all]

        csv (bool): input files are in tabular CSV format if True else input 
            is formatted as JSON. If input is in CSV format, field names will 
            be taken from the header (first line of the file).

    Returns:
        dict: relational database loaded as a dictionary
    """

    mapping = {}
    for map_file in infiles:
        in_h = open_io(map_file)

        if not csv:  #mapping is JSON formatted
            json_map = json.load(in_h)

            # Subset based on desired fields
            if fields:
                fields = set(fields)

                for item in json_map:
                    entry = json_map[item]
                    entry = {k: entry[k] for k in entry.keys() & fields}

                    # Add entry to DB
                    mapping[item] = entry
            else:
                # Add all entries to DB at once
                mapping = {**json_map, **mapping}

        else:  #mapping is CSV formatted
            header = in_h.readline()
            header = header.rstrip().split('\t')

            # Subset based on desired fields
            if not fields:
                fields = header

            # Extract indices corresponding to desired fields
            try:
                keep = [header.index(i) for i in fields]
            except ValueError:  #field not in header
                bad_fields = ', '.join(set(fields).difference(header))
                raise InputError("{}: line 1. Provided field(s) '{}' not "
                    "found in the file header".format(map_file, bad_fields))

            # Add entries to DB one per line
            for nline, line in enumerate(in_h):
                line = line.rstrip().split('\t')
                acc = line[0]
                try:
                    mapping[acc] = {header[i]:line[i] for i in keep}
                except IndexError:
                    raise FormatError("{}: line {}. The number of fields in "
                    "the row do not match the number of fields in the header"\
                    .format(map_file, nline))

    return(mapping)

def filter_dbs(mapping: dict, patterns: list, subset: bool=True, \
               case: bool=False):
    """Search for entries in the database matching one or more of the provided
    patterns

    Args:
        mapping (dict): dictionary containing the relational database entries
            to be filtered

        patterns (list): list of pattern matching criteria used to filter the 
            relational database. Each item should be of the form FIELD:PATTERN

        subset (bool): keep only those entries with values matching one or more
            of the provided patterns. If False, will only keep those entries 
            that do not match any of the provided patterns

        case (bool): case-sensitive pattern matching [default: False]

    Returns:
        dict: filtered relational database
    """

    # Speedup trick
    list_type = type(list())

    # Create dictionary for pattern-matching
    merged = {}
    for pattern in patterns:
        try:
            field, regex = pattern.split(':', 1)
        except ValueError:
            print("error: unable to parse search criteria {}".format(pattern), 
                  file=sys.stderr)
            sys.exit(1)

        if not case:
            regex = regex.lower()

        # Combine search patterns from the same field
        if field in merged:
            regex = "{}|{}".format(regex, merged[field])

        merged[field] = regex

    # Compile regular expressions
    criteria = {j: re.compile(k) for j, k in merged.items()}

    # Filter database using search criteria
    all_entries = list(mapping.keys())
    for entry_id in all_entries:
        match = False
        foi = set(criteria).intersection(mapping[entry_id])
        for field in foi:
            field_val = mapping[entry_id][field]
            if type(field_val) != list_type:
                field_val = [field_val]

            for val in field_val:
                if not case:
                    val = val.lower()

                match = criteria[field].search(val)
                if match:
                    break

            if match:
                break

        if not match and subset:
            del mapping[entry_id]
        elif match and not subset:
            del mapping[entry_id]

    return(mapping)

def derep_by_file(mapping: dict, inhandle):
    """Use input file to direct which entries should be combined into a 
    single entry

    Args:
        mapping (dict): dictionary containing the relational database entries
            to be dereplicated

        inhandle (filehandle): tab-separated file specifying which database 
            entries to combine

    Returns:
        dict: dereplicated relational database
    """
    for line in inhandle:
        if line.startswith('#'):
                continue

        split_line = line.strip().split('\t')

        # Combine entries from multiple databases
        merged = merge_entries(mapping, split_line)

        # Update database entries with combined information
        for entry in split_line:
            mapping[entry] = merged

    return(mapping)

def derep_by_field(mapping: dict, field: str):
    """Merge entries sharing the same field value
    Args:
        mapping (dict): dictionary containing the relational database entries
            to be dereplicated

        field (str): database field through which entries should be combined

    Returns:
        dict: dereplicated relational database
    """
    novalue = 0

    reverse_map = {}
    for entry_id in mapping:
        entry = mapping[entry_id]
        field_value = get_value_str(entry, field)

        if field_value == 'NA':
            novalue += 1
            continue

        try:
            reverse_map[field_value].append(entry_id)
        except KeyError:
            reverse_map[field_value] = [entry_id]

    for field_id in reverse_map:
        entry_ids = reverse_map[field_id]
        if len(entry_ids) > 1:  #found replicates
            merged = merge_entries(mapping, entry_ids)
        else:  #no replicates
            merged = mapping[entry_ids[0]]

        # Reduce memory usage by removing entries from mapping
        for entry_id in entry_ids:
            del(mapping[entry_id])

        del(merged[field])  #replicate field is now entry ID, so remove
        mapping[field_id] = merged
    
    if novalue:
        print("warning: there were {} entries with no value for field {}"\
            .format(novalue, field), file=sys.stderr)

    return(mapping)

def merge_entries(mapping: dict, entries: list):
    """Merge entries in a relational database
    Args:
        mapping (dict): dictionary containing the relational database

        entries (list): list of entry IDs contained within the database that 
            should be merged

    Returns:
        dict: new entry created by combining the fields of all entries in the 
            entries list
    """
    list_type = type(list())

    merged = {}
    for entry_id in entries:
        try:
            entry = mapping[entry_id]
        except KeyError:
            print("warning: entry '{}' not found in the combined "
                  "relational database".format(entry_id), file=sys.stderr)
            continue

        for field in entry:
            try:
                entry_value = entry[field]
            except KeyError:
                continue
            else:
                entry_type = type(entry_value)

            try:
                merged_value = merged[field]
            except KeyError:
                merged[field] = entry_value
                continue
            else:
                merged_type = type(merged_value)

            if entry_value and not merged_value:
                merged[field] = entry_value
                continue
            elif (merged_value and not entry_value) \
                or not (entry_value and merged_value):
                continue

            if entry_type == list_type and merged_type == list_type:
                merged_field = entry_value + merged_value
            elif entry_type == list_type and merged_type != list_type:
                entry_value.append(merged_value)
                merged_field = entry_value
            elif merged_type == list_type and entry_type != list_type:
                merged_value.append(entry_value)
                merged_field = merged_value
            else:
                if str(entry_value) in str(merged_value):
                    merged[field] = merged_value
                    continue
                elif str(merged_value) in str(entry_value):
                    merged[field] = entry_value
                    continue
                else:
                    merged_field = [entry_value, merged_value]

            merged[field] = list(set(merged_field))

    return merged

def entry_as_csv(ident: str, entry: dict, fields=None, sep='\t'):
    """Reformat database entry as CSV"

    Args:
        ident (str): entry ID

        entry (dict): dictionary entry from the relational database

        fields (list): list of fields to include in the output

    Returns:
        str: entry reformatted as CSV string
    """
    if not fields:
        add_fields = sorted(entry.keys())
    else:
        add_fields = fields

    field_values = [ident]
    for field in add_fields:
        field_val = get_value_str(entry, field)

        field_values.append(field_val)

    return sep.join(field_values)

def get_value_str(entry: dict, field: str):
    """Output an entries field value as a string

    Args:
        entry (dict): dictionary entry from the relational database

        field (str): field to extract value string from

    Returns:
        str: entry field value formatted as a string
    """
    list_type = type(list())
    str_type = type(str())

    # Escape reserved characters
    escape_map = {ord('"'): '%22',
                  ord(';'): '%3B',
                  ord('\t'): '%09',
                 }

    try:
        field_value = entry[field]
    except KeyError:  #feature or entry field not in database
        field_value = 'NA'
    except TypeError:  #database not formatted correctly
        print("error: dont know what to do with entry {} and field {}"\
            .format(entry, field), file=sys.stderr)
        sys.exit(1)

    value_type = type(field_value)
    if value_type == list_type:  # Attribute has multiple values
        field_value = ';'.join([i.translate(escape_map) for i in field_value])
    elif value_type == str_type:
        field_value = field_value.translate(escape_map)

    if not field_value:
        field_value = 'NA'

    return field_value
