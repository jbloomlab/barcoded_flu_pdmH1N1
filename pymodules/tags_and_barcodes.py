"""Python functions related to processing tags and barcodes."""


import collections

import pandas as pd


def extract_tags(readiterator,
                 cellbarcodes,
                 start,
                 end,
                 *,
                 cellbc_tag='CB',
                 umi_tag='UB',
                 min_frac=0.5,
                 ):
    """Get tag (or barcode) sequence for each cell barcode and UMI.
    
    Parameters
    -----------
    readiterator : iterator over pysam.AlignedSegment
        Read iterator of the type returned by `pysam.fetch`.
    cellbarcodes : set
        Valid cell barcodes, discard reads not in this set.
    start : int
        Start of tag in reference in Python 0-based indexing.
    end : int
        End of tag in reference.
    cellbc_tag : str
        The SAM tag for the cell barcode.
    umi_tag : str
        The SAM tag for the UMI.
    min_frac : float
        Only call tag if >= this frac of reads for that cell barcode and UMI
        agree on sequence. A value of 0.5 means that more than half the reads
        for the cell barcode / UMI must agree on the tag sequence.
        
    Returns
    -------
    pandas.DataFrame
        Columns are 'cell_barcode', 'UMI', and 'tag'. Entry in 'tag' is
        'ambiguous' if the tag does not meet the `min_frac` cutoff.
        
    """
    counts = collections.defaultdict(lambda: collections.defaultdict(
                    lambda: collections.defaultdict(int)))
    sites = set(range(start, end))  # tag sites
    nsites = len(sites)  # number of sites in tag
    
    # count all tag sequences for each cell barcode / UMI
    for read in readiterator:
        if not read.has_tag(cellbc_tag):
                continue  # no cell barcode
        cellbc = read.get_tag(cellbc_tag)
        if cellbc not in cellbarcodes:
            continue  # cell barcode not on whitelist
        if not read.has_tag(umi_tag):
            continue  # no UMI
        umi = read.get_tag(umi_tag)
        readseq = read.query_sequence
        tagseq = ''.join(readseq[iquery] for iquery, iref in
                         read.get_aligned_pairs(matches_only=True)
                         if iref in sites)
        if len(tagseq) != nsites:
            continue  # tag positions not fully covered
        counts[cellbc][umi][tagseq] += 1
        
    # now collapse to consensus tag sequences
    records = []
    for cellbc, umi_counts in counts.items():
        for umi, tag_counts in umi_counts.items():
            tot_counts = sum(tag_counts.values())
            assert tot_counts >= 1
            top_count, top_tag = sorted((count, tag) for tag, count in
                                        tag_counts.items())[-1]
            if top_count / tot_counts > min_frac:
                records.append((cellbc, umi, top_tag))
            else:
                records.append((cellbc, umi, 'ambiguous'))
    return pd.DataFrame.from_records(records,
                                     columns=['cell_barcode', 'UMI', 'tag'])
