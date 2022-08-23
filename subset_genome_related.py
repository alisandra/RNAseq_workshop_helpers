#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""make fasta/bam/gff files for sub-sequence of seq, with the range start_from-continue_to"""

from __future__ import print_function
### depends on samtools (tested 1.9)

from distutils.version import StrictVersion
import argparse
import os
import re
import subprocess
import random
import sys


# Named exceptions
class OutOfRangeError(Exception):
    pass


class DependencyIssuesError(Exception):
    pass


class BamWithoutFastaError(Exception):
    pass


# help parsing args
def parse_commas(x):
    """comma sep to python list for parsing user ARGS"""
    if x is None:
        out = []
    else:
        out = x.split(',')
    return out


# naming of output files
def get_name_out_fa(fasta, seq, start, end):
    """set output fasta name based on input name"""
    basename = re.sub('(\.fa$)|(\.fasta$)', '', fasta)
    return '{}__{}_{}-{}.fa'.format(basename, seq, start, end)


def get_name_out_bam(sam, seq, start, end):
    """set output bam name based on input name"""
    basename = re.sub('\.bam$', '', sam)
    return '{}__{}_{}-{}.bam'.format(basename, seq, start, end)


def get_name_out_gff(gff, seq, start, end):
    """set output gff name based on input name"""
    endings = re.match('.*(\.g[tf]f3?)$', gff)
    basename = re.sub(endings.group(1), '', gff)
    return '{}__{}_{}-{}{}'.format(basename, seq, start, end, endings.group(1))


# samtools versions are fun, but I don't ultimately know which all will work, so warnings
def check_samtools(try_anyways):
    """raises error if samtools absent or too old, unless user has set try_anyways"""
    # all the version checking necessary because samtools remains unstable in terms of parameters
    min_tested = "1.9" 
    max_tested = "1.9"
    if try_anyways:
        return 0
    samtools = 'samtools'
    version = None
    try:
        version_str = subprocess.check_output([samtools, '--version'])
    except OSError:
        raise DependencyIssuesError("Command {0} not recognized, {0} ({1}-{2}) must be installed/in $PATH".format(
            samtools, min_tested, max_tested
        ))

    for line in version_str.decode().split('\n'):
        match = re.match('.*samtools *([.0-9]*).*', line)
        if match:
            try:
                version = str(match.group(1))
            except ValueError:
                pass
    if version is None:
        raise DependencyIssuesError("""
{0} present but version could not be parsed from string
{1}
This script has been tested only for {0} {2}-{3}, if you wish to continue, please add --try_anyways""".format(
            samtools, version_str, min_tested, max_tested
        ))
    if StrictVersion(version) < StrictVersion(min_tested):
        
        raise DependencyIssuesError("""
{0} present but the identified version ({1}) is less than {2} where this script was tested
if you wish to continue, please add --try_anyways""".format(samtools, version, min_tested))
    elif StrictVersion(version) > StrictVersion(max_tested):
        print("WARN: {} newer ({}) than {} where this script was tested, should hopefully work, but be warned".format(
            samtools, version, max_tested
        ), file=sys.stderr)
    return 0


# actually parsing subsections
def get_range_fasta(fasta, seq, start, end):
    """creates and indexes fasta file of just the requested sub sequence (by calling samtools)"""
    fasta_out = get_name_out_fa(fasta, seq, start, end)
    print('cropping {} and writing to {}'.format(fasta, fasta_out))
    if not os.path.exists(fasta + '.fai'):
        subprocess.call(['samtools', 'faidx', fasta])
    # save subsequence to fasta_out
    subprocess.call(['samtools', 'faidx', fasta, '{}:{}-{}'.format(seq, start, end), '-o', fasta_out])
    subprocess.call(['samtools', 'faidx', fasta_out])  # and reindex
    return fasta_out


def get_range_gff(gff, seq, start, end):
    """creates gff file of just the requested subsequence and shifts coordinates"""
    gff_out = get_name_out_gff(gff, seq, start, end)
    print('cropping {} and writing to {}'.format(gff, gff_out))
    with open(gff_out, 'w') as fout:
        with open(gff) as fin:
            for line in fin:
                try:
                    fout.write(shift_gff_line(line, seq, start, end) + '\n')
                except OutOfRangeError:  # this is also the filter
                    pass


def get_range_bam(bam, seq, start, end, fai):
    """creates and indexes bam file of just the requested sub sequence and shifts coordinates"""
    bam_out = get_name_out_bam(bam, seq, start, end)
    print('cropping {} and writing to {}'.format(bam, bam_out))
    reads_flat = subprocess.check_output(['samtools', 'view', bam, '{}:{}-{}'.format(seq, start, end)])
    if len(reads_flat) == 0:
        reads = []
    else:
        reads = reads_flat.split(b'\n')
        if reads[-1] == '':
            reads.pop()
    tmp_sam = "%032x.sam" % random.getrandbits(128)
    with open(tmp_sam, 'w') as f:
        for line in reads:
            try:
                f.write(shift_sam_line(line, start, end) + '\n')
            except OutOfRangeError:  # throwing out partial overlap at start, otherwise crop sequence/cigar
                pass
    # todo, use exit stati
    subprocess.call(['samtools', 'view', '-ht', fai, '-b', tmp_sam, '-o', bam_out])
    subprocess.call(['samtools', 'index', bam_out])
    os.remove(tmp_sam)


def shifted_coordinates(start, stop, ori_start, ori_stop):
    """calculate local coordinates on sub-sequence start-stop for input global 'ori' coordinates"""
    # all param coordinates count from 1 and stop is included, bc bioinfo formats x_x
    shift_by = start - 1
    length = stop - shift_by
    out_start = ori_start - shift_by
    out_stop = ori_stop - shift_by
    if out_stop <= 0 or out_start > length:
        raise OutOfRangeError('{} {} out of new 1 - {} range'.format(out_start, out_stop, length))
    # the following only if some of range has overlap
    # trim ends
    out_start = max(1, out_start)
    out_stop = min(out_stop, length)
    return out_start, out_stop


def shift_gff_line(line, seq, start, stop):
    """adjust/shift one gff line from global coordinates to coordinates local to sub sequence"""
    line = line.rstrip()
    if line.startswith('#'):
        return line
    else:
        sline = line.split('\t')
    # update ref name
    if sline[0] != seq:
        raise OutOfRangeError('wrong chr {} not target {}'.format(sline[0], seq))
    sline[0] = '{}:{}-{}'.format(sline[0], start, stop)
    # update coordinates
    mod_start, mod_stop = shifted_coordinates(start, stop, ori_start=int(sline[3]), ori_stop=int(sline[4]))
    sline[3] = str(mod_start)
    sline[4] = str(mod_stop)

    return '\t'.join(sline)


def shift_sam_line(line, start, stop):
    """adjust/shift one sam line from global coordinates to coordinates local to sub sequence"""
    # note, this does not handle pairs falling out of range correctly, hopefully samtools view does
    line = line.rstrip()
    line = line.decode()
    if line.startswith('@'):
        return line
    elif not line:
        return line  # empty string
    else:
        sline = line.split('\t')

    # update 'ref' sequence to match samtools' output
    sline[2] = '{}:{}-{}'.format(sline[2], start, stop)

    try:
        ori_start = int(sline[3])
        length = int(sline[8])
    except ValueError:
        raise OutOfRangeError('read not mapped? {}'.format(line))

    # shift coordinates
    mod_start, mod_stop = shifted_coordinates(start, stop, ori_start=ori_start, ori_stop=ori_start + length)
    sline[3] = str(mod_start)  # pos
    try:
        pair_start = mod_start + (int(sline[7]) - ori_start)
        sline[7] = str(pair_start)
    except ValueError:
        pass  # e.g. without a pair there will be no number in sline[7]/PNEXT
    return '\t'.join(sline)


# and finally, flow control
def main(fasta, bams, gffs, seq, start_from, continue_to, try_anyways):
    """handles -> subsequence conversion of each provided file"""
    check_samtools(try_anyways)
    if fasta is None and bams is not None:
        raise BamWithoutFastaError("Cannot reconstruct bam header without fasta, please specify --fasta <your fasta>")

    if fasta is not None:
        fasta_out = get_range_fasta(fasta, seq, start_from, continue_to)

        for bam in parse_commas(bams):
            get_range_bam(bam, seq, start_from, continue_to, fasta_out + '.fai')

    for gff in parse_commas(gffs):
        get_range_gff(gff, seq, start_from, continue_to)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='tool to organize subsetting of fasta, bam, gff. All samtools sorting '
                                                 'and indexing of .bams must already be done. Ranges start at 1 and '
                                                 'are inclusive so as to match samtools / gff coordinates etc...')
    parser.add_argument('--fasta', nargs='?', default=None, help='fasta file to subset')
    parser.add_argument('--bam', nargs='?', default=None, help='bam file to subset (comma separate for multiple),'
                                                               'requires --fasta')
    parser.add_argument('--gff', nargs='?', default=None, help='gff file to subset (comma separate for multiple)')

    parser.add_argument('-s', '--seq', required=True, help='target sequence')
    parser.add_argument('-f', '--start', default=1, type=int, help='starting _from_ this bp (count from 1, because...)')
    parser.add_argument('-t', '--end', default=1e16, type=int, help='continue _to_ this bp')
    parser.add_argument('--try_anyways', action='store_true', help='ignores any errors/warnings on samtools versions')
    args = parser.parse_args()
    main(args.fasta, args.bam, args.gff, args.seq, args.start, args.end, args.try_anyways)
