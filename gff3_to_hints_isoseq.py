#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""converts gff3 of collapsed transcripts (AKA with just transcript and exon features) into hints for Augustus"""

from __future__ import print_function

from dustdas import gffhelper
import argparse
import sys


class Transcript(object):
    supported_targets = ('ass', 'dss', 'exon', 'intron', 'tss', 'tts', 'ep', 'exonpart', 'ip', 'intronpart')

    @property
    def known_targets(self):
        # if you add or remove targets here, they need to be added/removed for Transcript.supported_targets as well!
        out = {
            'ass': self.acceptor_splice_sites,
            'dss': self.donor_splice_sites,
            'exon': self.exons,
            'intron': self.introns,
            'tss': self.transcription_start,
            'tts': self.transcription_termination,
            'ep': self.exon_parts,
            'exonpart': self.exon_parts,
            'ip': self.intron_parts,
            'intronpart': self.intron_parts
        }

        return out

    def __init__(self, a_transcripts_worth, trim_exonparts=4, trim_intronparts=4):
        transcript = a_transcripts_worth.pop(0)
        # validate transcript set
        assert transcript.type == "transcript"
        assert transcript.start == a_transcripts_worth[0].start
        assert transcript.end == a_transcripts_worth[-1].end
        assert transcript.strand in ['-', '+']
        self.transcript = transcript
        self.group = self.transcript.get_ID()
        self.raw_exons = a_transcripts_worth

        for exon in self.raw_exons:
            if exon.type != "exon":
                raise ValueError("This script only knows how to handle 'transcript' and 'exon' features")
            assert exon.get_Parent()[0] == self.group
            exon.start = int(exon.start)
            exon.end = int(exon.end)

        self.exons = self._exons()
        self.introns = self._introns()
        self.transcription_start = self._transcription_start()
        self.transcription_termination = self._transcription_termination()
        self.acceptor_splice_sites = self._splice_sites(acceptor=True)
        self.donor_splice_sites = self._splice_sites(acceptor=False)
        self.exon_parts = self._parts(exons=True, trim=trim_exonparts)
        self.intron_parts = self._parts(exons=False, trim=trim_intronparts)

    def _exons(self):
        out = []
        for i in range(len(self.raw_exons)):
            out.append((self.raw_exons[i].start, self.raw_exons[i].end))
        return out

    def _introns(self):
        out = []
        for i in range(len(self.raw_exons) - 1):
            start_end = (self._lower_intron_end(i), self._higher_intron_end(i))
            out.append(start_end)
        return out

    def _splice_sites(self, acceptor=True):
        plus_strand = self.transcript.strand == "+"
        # take the 'higher number' when:
        if (plus_strand and acceptor) or (not plus_strand and not acceptor):
            site_fn = self._higher_intron_end
        else:
            site_fn = self._lower_intron_end
        out = []
        for i in range(len(self.raw_exons) - 1):
            at = site_fn(i)
            start_end = (at, at)
            out.append(start_end)
        return out

    def _higher_intron_end(self, i):
        return self.raw_exons[i + 1].start - 1

    def _lower_intron_end(self, i):
        return self.raw_exons[i].end + 1

    def _transcription_start(self):
        if self.transcript.strand == "+":
            at = self.transcript.start
        elif self.transcript.strand == "-":
            at = self.transcript.end
        else:
            raise ValueError("strand not in ['+', '-']")
        start_end = (at, at)
        return [start_end]  # keep it as a list of (start, end) tuples, so that all features have same type

    def _transcription_termination(self):
        if self.transcript.strand == "-":
            at = self.transcript.start
        elif self.transcript.strand == "+":
            at = self.transcript.end
        else:
            raise ValueError("strand not in ['+', '-']")
        start_end = (at, at)
        return [start_end]

    def _parts(self, exons=True, trim=4):
        if exons:
            wholes = self.exons
        else:
            wholes = self.introns
        out = []
        for start, end in wholes:
            try:
                out.append(self.trim_feature(start, end, trim))
            except CoordinateError:  # skips any 'parts' smaller than 2x trim
                pass
        return out

    @staticmethod
    def trim_feature(start, end, trim=4):
        """shrinks feature by 'trim' on both ends"""
        start += trim
        end -= trim

        if not start < end:
            raise CoordinateError
        return start, end


    def make_lines(self, targets, priority, source):
        """setup all hint gff3 lines for all requested features of a transcript"""
        known_targets = self.known_targets
        list_out = []
        for target in targets:
            feature_list = known_targets[target]

            for start, end in feature_list:
                list_out.append(self._make_line(start, end, target, priority=priority, source=source))
        return '\n'.join(list_out) + '\n'

    def _make_line(self, start, end, feature, priority=4, source="E"):
        """reformats into hints style gff3 line"""
        out = "{seqid}\t{this_script}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{phase}\tgrp={group};" \
              "pri={priority};src={source}".format(seqid=self.transcript.seqid,
                                                   this_script="gff3_to_hints_isoseq",
                                                   feature=feature,
                                                   start=start,
                                                   end=end,
                                                   score=self.transcript.score,
                                                   strand=self.transcript.strand,
                                                   phase=self.transcript.phase,
                                                   group=self.group,
                                                   priority=priority,
                                                   source=source)
        return out


def group_transcripts(gff):
    """generate group of gff lines corresponding to one transcript (so with [transcript, exon, exon, exon, ...])"""
    gh = gffhelper.read_gff_file(infile=gff)
    a_transcripts_worth = [next(gh)]
    for entry in gh:
        if entry.type == "transcript":
            yield a_transcripts_worth
            a_transcripts_worth = [entry]
        else:
            a_transcripts_worth.append(entry)
    yield a_transcripts_worth


class CoordinateError(Exception):
    pass


def validate_hint_types(user_hints):
    """parse and validate user input for hint types"""
    hints = user_hints.split(',')
    ok = True
    supported = Transcript.supported_targets
    for hint in hints:
        if hint not in supported:
            ok = False
            print('"{}" hint not recognized'.format(hint), file=sys.stderr)
    if not ok:
        raise ValueError("user provided hint(s) not recognized from: {}\n"
                         "acceptable values are: {}\nmultiple hints can be comma separated e.g. 'ass,dss,ip,ep'"
                         "".format(user_hints, supported))
    return hints


def main():
    # parse all the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--gff3_in', required=True, help='input gff3 file with "transcript" and "exon" features')
    parser.add_argument('-o', '--hints_out', required=True, help='output file')
    parser.add_argument('-p', '--priority', default=4, type=int,
                        help='priority with which Augustus should consider these hints, hints with higher priority will'
                             ' be preferred in case of conflicting hints by Augustus')
    parser.add_argument('-s', '--source', default='E', help='source info for Augustus')
    parser.add_argument('--trim_exonparts', default=4, type=int, help='trim exonpart (ep) hints by')
    parser.add_argument('--trim_intronparts', default=4, type=int, help='trim intronpart (ip) hints by')
    parser.add_argument('-t', '--hint_types', default='ass,dss,ep,ip,tss,tts',
                        help='comma separated list of hint types that should be produced. '
                             'Supported values are {}'.format(Transcript.supported_targets))
    args = parser.parse_args()

    hints = validate_hint_types(args.hint_types)

    with open(args.hints_out, 'w') as handleout:
        for a_transcripts_worth in group_transcripts(args.gff3_in):
            transcript = Transcript(a_transcripts_worth, trim_exonparts=args.trim_exonparts,
                                    trim_intronparts=args.trim_intronparts)
            handleout.write(transcript.make_lines(hints, priority=args.priority, source=args.source))


if __name__ == "__main__":
    main()
