"""Classify reads in a FASTQ file by termini and optionally trim.

Written by Jesse Bloom, 2019.

"""


import argparse
import contextlib
import math
import os

import matplotlib
import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import pysam

import regex


def parse_fastq_termini(fastq, fastqmate, multi_pattern, outprefix,
                        patterns, fastq_trim5, statsfile,
                        fastq_suffix='R1', fastqmate_suffix='R2',
                        no_match_name='no_match'):
    """Parse and classify reads by termini.

    Parameters
    ----------
    fastq : str
        FASTQ file with reads to parse and classify by termini.
    fastqmate : str
        FASTQ file with read mates of those in `fastq`.
    multi_pattern : {'first', 'combined'}
        If read termini matches multiple patterns, assign to first category
        or a combined category with all matches.
    outprefix : str
        Output prefix on created files. Directories created if needed.
    patterns : list
        A list of the 3-tuples `(name, pattern, mismatch)` giving pattern
        to match at beginning of `fastq` reads as a regex, the name of the
        pattern, and the maximum number of allowed substitution mismatches.
    fastq_trim5 : int
        Amount to trim from 5' end of reads in `fastq` that match a pattern.
    statsfile : str
        Name of created CSV file with parsing stats.
    fastq_suffix : str
        Suffix on parsed `fastq` output.
    fastqmate_suffix : str
        Suffix on parsed `fastqmate` output.
    no_match_name : str
        Name for output file for reads that match no pattern in `patterns`.

    Note
    ----
    The created output files have names of the form:
    ``{outprefix}_{name}_{suffix}.fastq`` where ``{name}`` is the name of
    the pattern that matches, and ``{suffix}`` is `fastq_suffix` or
    `fastqmate_suffix`.

    """

    pattern_name = [tup[0] for tup in patterns]
    pattern = [tup[1] for tup in patterns]
    pattern_mismatch = [tup[2] for tup in patterns]

    if len(pattern_name) != len(set(pattern_name)):
        raise ValueError('duplicate pattern names')
    if no_match_name in pattern_name:
        raise ValueError(f"`no_match_name` {no_match_name} is a pattern name ")
    if fastq_suffix == fastqmate_suffix:
        raise ValueError('`fastq_suffix` and `fastqmate_suffix` are the same')
    if fastq_trim5 < 0:
        raise ValueError('`fastq_trim5` must be >= 0')
    if multi_pattern not in {'first', 'combined'}:
        raise ValueError(f"invalid `multi_pattern` {multi_pattern}")

    # compile regex with fuzzy matching allowing substitutions:
    # https://pypi.org/project/regex/
    pattern_regexs = {name: regex.compile(f"(?:{pattern}){{s<={mismatch}}}")
                      for name, pattern, mismatch in zip(pattern_name,
                                                         pattern,
                                                         pattern_mismatch)}

    if os.path.dirname(outprefix):
        os.makedirs(os.path.dirname(outprefix), exist_ok=True)

    outfilenames = []
    outfiles = {}
    counts = {}
    matches_name = {}
    with contextlib.ExitStack() as stack:
        # Define callback to delete files on error. See here:
        # https://docs.python.org/3/library/contextlib.html#replacing-any-use-of-try-finally-and-flag-variables
        @stack.callback
        def delete_files_on_err():
            for fname in outfilenames:
                if os.path.isfile(fname):
                    os.remove(fname)

        # open FASTQ files to read
        fastq = stack.enter_context(pysam.FastxFile(fastq))
        fastqmate = stack.enter_context(pysam.FastxFile(fastqmate))

        # iterate over all read pairs
        for r, rmate in zip(fastq, fastqmate):

            if r.name != rmate.name:
                raise IOError(f"read name mismatch:\n{r.name}\n{rmate.name}")

            seq = r.sequence

            # get patterns matched by this read
            matches = []
            for name, pattern_regex in pattern_regexs.items():
                if pattern_regex.match(seq):
                    matches.append(name)
                    if multi_pattern == 'first':
                        break
            matches = tuple(matches)

            # if first read matching this pattern, set up counts and files
            if matches not in counts:
                counts[matches] = 0
                if not matches:
                    matches_name[matches] = no_match_name
                else:
                    matches_name[matches] = '_and_'.join(matches)
                for f, suffix in [('fastq', fastq_suffix),
                                  ('fastqmate', fastqmate_suffix)]:
                    fname = f"{outprefix}_{matches_name[matches]}_{suffix}.fastq"
                    outfilenames.append(fname)
                    outfiles[(matches, f)] = stack.enter_context(open(fname, 'w'))

            # add to counts and write to appropriate outfile
            counts[matches] += 1
            if matches:
                # trim if a match
                seq = seq[fastq_trim5:]
                qual = r.quality[fastq_trim5:]
            else:
                qual = r.quality
            outfiles[(matches, 'fastq')].write(
                    f"@{r.name} {r.comment}\n{seq}\n+\n{qual}\n")
            outfiles[(matches, 'fastqmate')].write(f"{rmate}\n")

        # make sure both FASTQ files at end of file
        for f in [fastq, fastqmate]:
            try:
                next(fastq)
            except StopIteration:
                pass
            else:
                raise IOError(f"Not at end of both `fastq` {fastq} and "
                              f"`fastqmate` {fastqmate}. Maybe these "
                              'files have different numbers of entries?')

        # make CSV file with stats
        counts = {matches_name[matches]: n for matches, n in counts.items()}
        if statsfile in outfiles:
            raise ValueError(f"illegal `statsfile` name {statsfile}")
        (pd.Series(counts)
         .rename('n_reads')
         .rename_axis('pattern')
         .sort_values(ascending=False)
         .to_csv(statsfile, header=True)
         )

        stack.pop_all()  # no callback to delete files if reach here


def plot_fastq_termini_stats(statsfiles, plotfile, run_names, aggregate,
                             maxcol=6):
    """Plot summary of stats from :func:`parse_fastq_termini`.

    Parameters
    ----------
    statsfiles : list
        Stats files created by :func:`parse_fastq_termini`.
    plotfile : str
        Name of created plot file.
    run_names : list or `None`.
        Run names corresponding to each file in `statsfiles`. Can be `None`
        if `aggregate` is `True`.
    aggregate : bool
        Aggregate statistics across all runs, or facet plot by run.
    maxcol : int
        Max number of columns in faceted plot.

    """
    if aggregate:
        stats = pd.concat(pd.read_csv(f) for f in statsfiles)
        ncol = nrow = 1
    else:
        if len(statsfiles) != len(run_names):
            raise ValueError('`statsfiles` and `run_names` differ in length')
        stats = pd.concat(pd.read_csv(f).assign(run=run) for f, run in
                          zip(statsfiles, run_names))
        ncol = min(len(statsfiles), maxcol)
        nrow = math.ceil(len(statsfiles) / ncol)

    p = (ggplot(stats,
                aes('pattern', 'n_reads', label='n_reads', fill='pattern')) +
         geom_bar(stat='identity') +
         geom_text(va='bottom', size=7) +
         theme(axis_text_x=element_text(angle=90),
               figure_size=(3 * ncol, 3 * nrow),
               legend_position='none',
               )
         )

    if not aggregate:
        p = p + facet_wrap('~ run', ncol=ncol)

    backend = matplotlib.get_backend()
    plt.switch_backend('Cairo')
    p.save(plotfile, verbose=False)
    plt.switch_backend(backend)
