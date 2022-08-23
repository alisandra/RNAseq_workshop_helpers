#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Find and print all AGIs (via RegEx) in a text file"""

from __future__ import print_function

import re
import sys
import getopt


def get_all_agis(file_in, re_agi):
    """Find all AGIs in a text and ignore all context"""
    text_in = open(file_in).read()
    matches_agi = re_agi.findall(text_in)
    matches_agi = [refomat_agi(x) for x in matches_agi]
    return matches_agi


def get_agis_by_line(file_in, re_agi, delimiter='\t', sub_delimiter=';'):
    """Find all AGIs per line and report these before each line"""
    lines_in = open(file_in).readlines()
    out = []
    for line in lines_in:
        line = line.rstrip()
        matches_agi = re_agi.findall(line)
        matches_agi = [refomat_agi(x) for x in matches_agi]
        matches_str = sub_delimiter.join(matches_agi)
        newline = matches_str + delimiter + line
        out.append(newline)
    return out


def refomat_agi(agi):
    out = agi.upper()
    return out


def usage():
    usagestr = """ python agi_finder.py -i text_file [options] > AGIs.txt
###############
-i | --in=              input text file
-l | --line_wise        find AGI #s by line and return with line
-h | --help             prints this message
"""
    print(usagestr, file=sys.stderr)
    sys.exit(1)


def main():
    file_in = None
    by_line = False
    # get opt
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:lh", ["in=", "line_wise", "help"])
    except getopt.GetoptError as err:
        print (str(err), file=sys.stderr)
        usage()

    for o, a in opts:
        if o in ("-i", "--in"):
            file_in = a
        elif o in ("-l", "--line_wise"):
            by_line = True
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if file_in is None:
        print("input file required")
        usage()

    # RegEx by which AGI's are actually identified
    re_agi = re.compile('[Aa][Tt][CcMm1-5][Gg][0-9]{5}')

    if by_line:
        to_print = get_agis_by_line(file_in, re_agi)
    else:
        to_print = get_all_agis(file_in, re_agi)

    for unit in to_print:
        print(unit)


if __name__ == "__main__":
    main()
