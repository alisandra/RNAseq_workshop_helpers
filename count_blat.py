#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""count and report hits to each target in user provided blat output (.psl)"""
from __future__ import print_function

import getopt
import sys


def usage(msg=''):
    usestr = """count_blat.py -i blat_file.psl [options] > blat_counts.tsv
counts the hits from a blast/blat file that match each target
###########################################################
requires:
-i|--in=        a blast (m8/outfmt6) or blat (out=blast8) formatted file

optional:
-s|--sorted     if the file is already sorted by query and score this reduces required computation/memory
-a|--all        set this parameter if you wish to count -all- and not just -best- blat hits
-h|--help       prints this
"""
    print(usestr + '\n' + msg)
    sys.exit(1)


class BlastCounter:
    q_id_col = 0
    t_id_col = 1
    score_col = 11

    def by_query(self, iterable):
        """generate chunks of the blat output with matching query ID"""
        prev_query = None
        by_query = []
        for line in iterable:
            line = line.rstrip()
            line = line.split()
            query = line[self.q_id_col]
            target = line[self.t_id_col]
            score = float(line[self.score_col])
            if query == prev_query or prev_query is None:
                by_query.append([query, target, score])
            else:
                yield by_query
                by_query = [[query, target, score]]
            prev_query = query

    def count(self, filein, re_sort=True, best_only=True):
        """count how many hits there are to each target sequence in blat output"""
        seen = {}
        by_targets = {}

        # sorting by query name (first thing in each line)
        if re_sort:
            lines = open(filein).readlines()
            lines = sorted(lines)
        # if pre-sorted
        else:
            lines = open(filein)
        for query_set in self.by_query(lines):
            # confirm sort by query ID
            if query_set[0][0] in seen:
                raise SortingError("{} has been seen in non-consecutive blocks. \n\
                If -s is set, input must be presorted by query ID.".format(query_set[0][0]))
            else:
                seen[query_set[0][0]] = 1
            # sort by bit score if necessary
            if re_sort:
                query_set = sorted(query_set, key=lambda a_hit: a_hit[2], reverse=True)
            # keep only top hits
            if best_only:
                query_set = [query_set[0]]
            # add remaining hits to their target sequences
            for hit in query_set:
                self.init_or_incr(by_targets, hit[1])

        if isinstance(lines, file):
            lines.close()

        return by_targets

    @staticmethod
    def init_or_incr(dictionary, key):
        if key in dictionary:
            dictionary[key] += 1
        else:
            dictionary[key] = 1


class SortingError(Exception):
    pass


def main():
    """interpret user input, count and report hits to each target sequence"""
    # default parameters
    best = True
    re_sort = True
    filein = None

    # this whole section interprets the command line parameters
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "i:sah",
                                       ["in=", "sorted", "all", "help"])
    except getopt.GetoptError as err:
        print (str(err))
        usage()

    for o, a in opts:
        if o in ("-i", "--in"):
            filein = a
        elif o in ("-s", "--sorted"):
            re_sort = False
        elif o in ("-a", "--all"):
            best = False
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if filein is None:
        usage("input blast file required (-i)")

    # process the blat file
    # all the mechanics are found in the class BlatCounter
    blast_counter = BlastCounter()
    counts_by_target = blast_counter.count(filein, re_sort=re_sort, best_only=best)
    if not counts_by_target:
        print('WARN: no hits found, is the input file empty?', file=sys.stderr)
    # counts_by_target is now a dictionary with target IDs as keys and number of hits as values
    for target in counts_by_target:
        print ('{}\t{}'.format(target, counts_by_target[target]))


if __name__ == "__main__":
    main()
