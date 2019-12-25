"""Utilities for plotting information associated with viral genes."""


import itertools
import re

import Bio.Data.IUPACData

import dna_features_viewer

import matplotlib
import matplotlib.pyplot as plt

import pandas as pd


def plot_genes_and_coverage(genes,
                            coverage_df,
                            *,
                            genes_relheight=0.3,
                            gene_features={'exon': '#999999', 
                                           'viral_barcode': '#E69F00',
                                           'viral_tag': '#56B4E9',
                                           },
                            full_gene_feature='exon',
                            site_tick_interval=250,
                            color_mutations=False,
                            mutation_colors={'A': '#009E73',
                                             'C': '#F0E442',
                                             'G': '#0072B2',
                                             'T': '#D55E00',
                                             },
                            figwidth=17,
                            ):
    """Plot gene structure and coverage for one or more samples.
    
    Parameters
    -----------
    genes : list
        `BioPython SeqRecord <https://biopython.org/wiki/SeqRecord>`_ for
        each gene.
    coverage_df : pandas.DataFrame
        Coverage at each site. Must have columns 'gene', 'sample', 'site', and
        'coverage'. Entries in 'gene' column should match `id` attributes
        of the SeqRecords in `genes`. Sites should be in 1, ... numbering.
    genes_relheight : float
        Relative height of top row with genes relative to coverage plots.
    gene_features : dict
        Which features in `genes` to plot. Keyed by regex strings that match
        `type` attributes of SeqRecords in `genes`; values are colors to plot.
        Features with types not in `gene_features` are not plotted.
    full_gene_feature : str
        Name of feature in each gene for full gene length. Must be in
        `gene_features` too.
    site_tick_interval : int
        Frequency of ticks on x-axis giving site numbers.
    color_mutations : bool
        If `True`, overlay plots showing counts of each mutation at each site.
        Mutation colors are **not** stacked, but overlaid on the coverage and
        each other. Requires `coverage_df` to have columns giving reads at
        each site for each nucleotide. Only colors **mutant** identities;
        nucleotide that match the wildtype identity in `genes` (which can
        be ambiguous) are not colored.
    mutation_colors : dict
        Colors for mutations if using `color_mutations`.
    figwidth : float
        Overall figure width in inches.
    
    Returns
    --------
    matplotlib.figure.Figure, matplotlib.axes.Axes
        The figure and axes grid.
    
    """
    # check validity of genes and coverage
    gene_names = {gene.id for gene in genes}
    if len(genes) != len(set(gene_names)):
        raise ValueError('genes in `genes` do not have unique IDs')
    if set(coverage_df['gene']) != set(gene_names):
        raise ValueError('`coverage_df` and `genes` have different genes')
    
    # set up grid of subplots
    gene_lengths = [len(gene) for gene in genes]
    samples = list(coverage_df['sample'].unique())
    fig, axes = plt.subplots(
                    nrows=1 + len(samples),
                    ncols=len(genes),
                    sharex='col',
                    gridspec_kw={'width_ratios': gene_lengths,
                                 'height_ratios': ([genes_relheight] +
                                                   [1] * len(samples))
                                 },
                    figsize=(figwidth,
                             1.5 * (len(samples) + genes_relheight) + 1.5),
                    squeeze=False,
                    )

    # Convert genes to plottable graphic records
    class ViralGeneTranslator(dna_features_viewer.BiopythonTranslator):
        """Translate BioPython SeqRecord into GraphicRecord.""" 
        def compute_feature_color(self, feature):
            for regex, color in gene_features.items():
                if re.search(regex, feature.type):
                    return color

        def compute_feature_label(self, feature):
            return None
    
        def compute_filtered_features(self, features):
            return [f for f in features if
                    any(re.search(regex, f.type)
                        for regex in gene_features.keys())
                    ]
    
    # plot the genes in the first row of the subplots
    for ax, seqrecord in zip(axes[: 1,].ravel(), genes):
        graphic_record = ViralGeneTranslator().translate_record(seqrecord)
        for f in graphic_record.features:
            f.linecolor = f.color
        graphic_record.plot(ax=ax, with_ruler=False, draw_line=False)
        ax.set_title(seqrecord.id, fontsize=15)

    # plot coverage on subplots not in first row
    full_gene_feature_regex = [regex for regex in gene_features.keys()
                               if re.search(regex, full_gene_feature)]
    if len(full_gene_feature_regex) != 1:
        raise ValueError('no `full_gene_feature` in `gene_features`')
    else:
        full_gene_feature_regex = full_gene_feature_regex[0]
    for ax, (sample, seqrecord) in zip(axes[1:,].ravel(),
                                  itertools.product(samples, genes)):
        gene = seqrecord.id
        
        # data frame of coverage for each site by nt, whether it is wildtype
        df = (coverage_df
              .query('(sample == @sample) & (gene == @gene)')
              .melt(id_vars=['sample', 'gene', 'site'],
                    var_name='nucleotide',
                    value_name='coverage')
              .merge(pd.DataFrame.from_records(
                      enumerate(str(seqrecord.seq), start=1),
                      columns=['site', 'wildtype']),
                     on='site')
              .assign(
                  wildtype=lambda x: 
                           (x['wildtype']
                            .map(Bio.Data.IUPACData.ambiguous_dna_values)
                            ),
                  is_wildtype=lambda x: x.apply(lambda r: (r['nucleotide']
                                                           in r['wildtype']),
                                                axis=1)
                  )
              )
        
        # plot coverage over full gene
        (df
         .query('nucleotide == "coverage"')
         .plot(ax=ax,
               x='site',
               y='coverage',
               kind='area',
               legend=False,
               color=gene_features[full_gene_feature_regex],
               )
         )
        
        # optionally overlay a plot for mutations
        if color_mutations:
            (df
             .query('nucleotide != "coverage"')
             .query('not is_wildtype')
             .pivot_table(index='site',
                          columns='nucleotide',
                          values='coverage')
             .fillna(0)
             [list(mutation_colors)]
             .reset_index()
             .plot(ax=ax,
                   x='site',
                   kind='area',
                   stacked=False,
                   color=list(mutation_colors.values()),
                   #lw=1,
                   legend=False,
                   )
             )
            
        # set axes formatting
        ax.set_ylabel(None)
        ax.set_xlabel(None)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0),
                            useMathText=True)
        if sample == samples[-1]:
            xticklocs = list(range(1, int(ax.get_xlim()[1]),
                                   site_tick_interval))
            ax.set_xticks(ticks=xticklocs)
            ax.tick_params(axis='x', labelrotation=90)
        else:
            ax.set_xticklabels([])

        # plot rectangles for features of interest
        ymin, ymax = ax.get_ylim()
        for f in seqrecord.features:
            for regex, color in gene_features.items():
                if re.search(regex, f.type) and (regex !=
                                                 full_gene_feature_regex):
                    break
            else:
                continue
            start = f.location.start
            end = f.location.end
            rect = matplotlib.patches.Rectangle(
                                xy=(start - 0.5, ymin),
                                width=end - start + 1,
                                height=ymax,
                                color=color,
                                zorder=4,
                                alpha=0.4)
            ax.add_patch(rect)
        
    # shared x-, y-labels following here: https://stackoverflow.com/a/53172335
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False,
                    right=False)
    plt.xlabel('\nsite', size=14)
    plt.ylabel('coverage          ', size=14)
    
    # add legend for features
    feature_legend = plt.legend(
              handles=[matplotlib.patches.Patch(facecolor=color,
                                                 edgecolor=color,
                                                 label=regex)
                        for regex, color in gene_features.items()
                        if regex != full_gene_feature_regex],
               bbox_to_anchor=(1.006, ((len(samples) + 0.6 * genes_relheight) /
                                       (len(samples) + genes_relheight))),
               bbox_transform=plt.gcf().transFigure,
               fontsize=13,
               handlelength=0.5,
               title='features',
               title_fontsize=14,
               )
    
    # add legend for mutation colors
    if color_mutations:
        plt.gca().add_artist(feature_legend)
        plt.legend(handles=[matplotlib.patches.Patch(facecolor=color,
                                                     edgecolor=color,
                                                     label=nt)
                            for nt, color in mutation_colors.items()],
               bbox_to_anchor=(1, ((len(samples) - 0.2) /
                                       (len(samples) + genes_relheight))),
               bbox_transform=plt.gcf().transFigure,
               fontsize=13,
               handlelength=0.5,
               title='mutation',
               title_fontsize=14,
               ncol=2,
               )
    
    fig.tight_layout(w_pad=0)
    return fig, axes
