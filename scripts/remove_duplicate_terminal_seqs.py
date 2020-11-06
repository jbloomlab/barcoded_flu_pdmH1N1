"""Remove one copy of duplicated sequences termini of FASTQ reads..

Written by Jesse Bloom, 2020.

"""


import argparse
import gzip
import os

import pandas as pd

import pysam

import regex


def main():
    """Main body of script."""
    parser = argparse.ArgumentParser(
                description='Remove one copy of duplicated terminal seqs.',
                )
    parser.add_argument('fastq_in', help='Input FASTQ (can be gzipped).')
    parser.add_argument('fastq_out', help='Output FASTQ (can be gzipped).')
    parser.add_argument('terminal_seqs',
                        help='TSV or CSV with terminal sequences as second '
                             'column. Removes one of each duplicated copy '
                             'of any of these at either termini.')
    parser.add_argument('--n_mismatches', default=0,
                        help='Allow this many mismatches or indels in each '
                             'repeat of terminal sequences.')
    args = parser.parse_args()

    print(f"Reading terminal sequences from {args.terminal_seqs}")
    if os.path.splitext(args.terminal_seqs)[1] == '.tsv':
        terminal_seqs = pd.read_csv(args.terminal_seqs, sep='\t')
    else:
        terminal_seqs = pd.read_csv(args.terminal_seqs)
    if len(terminal_seqs.columns) != 2:
        raise ValueError(f"not exactly two columns in {args.terminal_seqs}")
    terminal_seqs = terminal_seqs[terminal_seqs.columns[1]].tolist()
    print(f"Read the following {len(terminal_seqs)} terminal sequences:\n  " +
          '\n  '.join(terminal_seqs))

    n = args.n_mismatches
    print(f"Allowing {n} single nt mismatches or indels.")
    match_str = []
    for terminal in terminal_seqs:
        match_str.append(f"({terminal}){{e<={n}}}({terminal}){{e<={n}}}")
    match_str = '(?:' + '|'.join(match_str) + ')'
    match5 = regex.compile('^' + match_str)
    match3 = regex.compile(match_str + '$')

    ntot = n5 = n3 = 0
    print(f"\nParsing reads in {args.fastq_in} to remove most terminal "
          "of duplicated termini, writing to {args.fastq_out}...")
    with pysam.FastxFile(args.fastq_in) as fastq_in, \
         gzip.open(args.fastq_out, mode='wt') as fastq_out:
        for read in fastq_in:
            ntot += 1
            m5 = match5.search(read.sequence)
            if m5:
                n5 += 1
                for first_termini in m5.groups():
                    if first_termini is not None:
                        read.sequence = read.sequence[len(first_termini):]
                        read.quality = read.quality[len(first_termini):]
                        break
                else:
                    raise RuntimeError('no matching group')
            m3 = match3.search(read.sequence)
            if m3:
                n3 += 1
                for last_termini in reversed(m3.groups()):
                    if last_termini is not None:
                        read.sequence = read.sequence[: -len(last_termini)]
                        read.quality = read.quality[: -len(last_termini)]
                        break
                else:
                    raise RuntimeError('no matching group')
            fastq_out.write(str(read))
            fastq_out.write('\n')
    print(f"Parsed {ntot} reads, removed {n5} and {n3} duplicate 5' and 3' "
          'termini, respectively.')


if __name__ == '__main__':
    main()
