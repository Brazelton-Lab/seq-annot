#! /usr/bin/env python
"""Functions for manipulating the relational databases
"""

import json
import re
from seq_annot.seqio import open_io, write_io
import sys


def load_dbs(infiles:list, fields: list=None):
    """Load relational databases into memory

    Args:
        infiles (list): list of input relational databases in JSON format

        fields (list): list of fields to load [default: load all]

    Returns:
        dict: relational database loaded as a dictionary
    """

    mapping = {}
    for map_file in infiles:
        json_map = json.load(open_io(map_file))
        if fields:
            for item in json_map:
                entry = json_map[item]
                entry = {k: entry[k] for k in entry.keys() & set(fields)}

                mapping[item] = entry
        else:
            mapping = {**json_map, **mapping}

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
            print("error: unable to parse search criteria", file=sys.stderr)
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
    reverse_map = {}
    for entry in mapping:
        try:
            rep = mapping[entry][field]
        except KeyError:
            print("error: field '{}' not found in the combined relational "
                  "database for entry '{}'".format(field, entry), \
                  file=sys.stderr)
            sys.exit(1)

        try:
            reverse_map[rep].append(entry)
        except KeyError:
            reverse_map[rep] = [entry]

    for rep in reverse_map:
        entries = reverse_map[rep]
        if len(entries) > 1:  #found replicates
            merged = merge_entries(mapping, entries)
        else:  #no replicates
            merged = mapping[entries[0]]

        # Reduce memory usage by removing entries from mapping
        for entry in entries:
            del(mapping[entry])

        del(merged[field])  #replicate field is now entry ID, so remove
        mapping[rep] = merged

    return(mapping)

def merge_entries(mapping: dict, entries: list):
    """Merge entries in a relational database
    Args:
        mapping (dict): dictionary containing the relational database

        entries (list): list of entries in database to combine together

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
