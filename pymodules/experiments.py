"""Modules for parsing experiments configuration."""


import copy

import pandas as pd

import ruamel.yaml


class Experiments:
    """Parse and access experiments configuration.

    Parameters
    ----------
    experiments_config : str or dict
        Name of YAML file with experiment configuration, or dict
        read from such a file.
    viral_tags : list
        Valid viral tags, parsed from yaml file specified in data
        directory.
    barcoded_viral_genes : list
        Viral genes with barcode, parsed from viral genbank file.

    Attributes
    ----------
    config_dict : dict
        Direct result of parsing the map in `experiments_config_file`.
    experiments : list
        All experiments.
    transcriptomic_runs : list
        All transcriptomic runs.
    transcriptomics_df : pandas.DataFrame
        Data frame of all transcriptomics runs.
    pacbio_df : pandas.DataFrame
        Data frame of all pacbio runs.
    expts_with_pacbio : list
        Experiments that have PacBio viral sequencing.
    viral_barcodes_df : pandas.DataFrame
        Data frame with all viral barcode sequencing paths
    expts_with_progeny_barcodes: list
        Experiments that have progeny viral barcode data

    """

    def __init__(self, experiments_config, viral_tags, barcoded_viral_genes):
        """See main class docstring."""
        if isinstance(experiments_config, dict):
            self.config_dict = copy.deepcopy(experiments_config)
        else:
            yaml = ruamel.yaml.YAML(typ='safe')
            with open(experiments_config) as f:
                self.config_dict = yaml.load(f)
            assert isinstance(self.config_dict, dict)

        self.experiments = list(self.config_dict)
        self._expect_ncells = {}

        # now get data for each experiment
        valid_keys = {'description',
                      'lab_notes',
                      'expect_ncells',
                      'transcriptomics',
                      'viral_barcodes',
                      'pacbio_viral_sequencing'}
        transcriptomics_records = []
        viral_barcodes_records = []
        pacbio_records = []
        for expt, expt_d in self.config_dict.items():

            if not set(expt_d).issubset(valid_keys):
                raise ValueError(f"invalid entries for {expt}")

            if 'transcriptomics' in expt_d:
                for run, run_d in expt_d['transcriptomics'].items():
                    transcriptomics_records.append((
                            expt,
                            run,
                            f"{expt}_{run}",
                            run_d['index'],
                            run_d['bcl_folder'],
                            run_d['lane'],
                            run_d['index_sequencing'],
                            ))

            if 'pacbio_viral_sequencing' in expt_d:
                for run, subread_bam in (expt_d['pacbio_viral_sequencing']
                                         .items()):
                    pacbio_records.append((
                            expt,
                            run,
                            f"{expt}_{run}",
                            subread_bam,
                            ))

            if 'expect_ncells' in expt_d:
                self._expect_ncells[expt] = expt_d['expect_ncells']
                if not isinstance(self._expect_ncells[expt], int):
                    raise ValueError(f"`expect_ncells` not int for {expt}")
            else:
                raise KeyError(f"`expect_ncells` missing for {expt}")

            if 'viral_barcodes' in expt_d:
                for source, source_d in expt_d['viral_barcodes'].items():
                    for tag, tag_d in source_d.items():
                        if tag not in viral_tags:
                            raise ValueError(f"invalid tag entry for {tag}")
                        for gene, gene_d in tag_d.items():
                            if gene not in barcoded_viral_genes:
                                raise ValueError(f"gene {gene} not barcoded")
                            for replicate, replicate_d in gene_d.items():
                                for run, fastq_path in replicate_d.items():
                                    viral_barcodes_records.append((
                                                                expt,
                                                                source,
                                                                tag,
                                                                gene,
                                                                replicate,
                                                                run,
                                                                fastq_path))

        self.transcriptomics_df = pd.DataFrame(transcriptomics_records,
                                               columns=['experiment',
                                                        'run_name',
                                                        'transcriptomic_run',
                                                        'index',
                                                        'bcl_folder',
                                                        'lane',
                                                        'index_sequencing',
                                                        ])
        self.transcriptomic_runs = (self.transcriptomics_df
                                    ['transcriptomic_run']
                                    .tolist()
                                    )
        assert (len(self.transcriptomic_runs) ==
                len(set(self.transcriptomic_runs)))

        self.pacbio_df = pd.DataFrame(pacbio_records,
                                      columns=['experiment',
                                               'run_name',
                                               'pacbio_run',
                                               'subread_bam',
                                               ])
        self.expts_with_pacbio = (self.pacbio_df
                                  ['experiment']
                                  .unique()
                                  .tolist()
                                  )
        assert len(self.pacbio_df) == self.pacbio_df['pacbio_run'].nunique()

        self.viral_barcodes_df = pd.DataFrame(viral_barcodes_records,
                                              columns=['experiment',
                                                       'source',
                                                       'tag',
                                                       'gene',
                                                       'replicate',
                                                       'run',
                                                       'fastq_path'])
        self.expts_with_progeny_barcodes = (self.viral_barcodes_df
                                            ['experiment']
                                            .unique()
                                            .tolist()
                                            )
        assert (len(self.expts_with_progeny_barcodes) ==
                len(set(self.expts_with_progeny_barcodes)))

    def transcriptomic_index(self, transcriptomic_run):
        """str: Illumina index for `transcriptomic_run`."""
        assert transcriptomic_run in self.transcriptomic_runs, (
                f"invalid `transcriptomic_run` {transcriptomic_run}")
        return str(self.transcriptomics_df
                   .set_index('transcriptomic_run')
                   .at[transcriptomic_run, 'index']
                   )

    def transcriptomic_bcl_folder(self, transcriptomic_run):
        """str: BCL folder for `transcriptomic_run`."""
        assert transcriptomic_run in self.transcriptomic_runs, (
                f"invalid `transcriptomic_run` {transcriptomic_run}")
        return str(self.transcriptomics_df
                   .set_index('transcriptomic_run')
                   .at[transcriptomic_run, 'bcl_folder']
                   )

    def transcriptomic_lane(self, transcriptomic_run):
        """str: Lane for `transcriptomic_run`."""
        assert transcriptomic_run in self.transcriptomic_runs, (
                f"invalid `transcriptomic_run` {transcriptomic_run}")
        return str(self.transcriptomics_df
                   .set_index('transcriptomic_run')
                   .at[transcriptomic_run, 'lane']
                   )

    def transcriptomic_index_sequencing(self, transcriptomic_run):
        """str: Index seq strategy for `transcriptomic_run`."""
        assert transcriptomic_run in self.transcriptomic_runs, (
                f"invalid `transcriptomic_run` {transcriptomic_run}")
        return str(self.transcriptomics_df
                   .set_index('transcriptomic_run')
                   .at[transcriptomic_run, 'index_sequencing']
                   )

    def expt_transcriptomic_runs(self, expt):
        """list: Transcriptomic runs for experiment `expt`."""
        assert expt in self.experiments, f"invalid `expt` {expt}"
        return (self.transcriptomics_df
                .query('experiment == @expt')
                ['transcriptomic_run']
                .tolist()
                )

    def expect_ncells(self, expt):
        """int: Expected number of cells for experiment `expt`."""
        return self._expect_ncells[expt]

    def expt_viral_barcode_fastqs(self, expt):
        """pandas.DataFrame: FASTQs of viral barcodes

        This function takes an argument with the desired
        experiment, `expt`. It returns a pandas
        DataFrame with paths to the viral barcode FASTQ
        files for that experiment. Each FASTQ path is
        annotated with source, tag, gene, replicate, and run.
        """
        assert expt in self.experiments, f"invalid `expt` {expt}"
        return (self.viral_barcodes_df
                .query('experiment == @expt')
                )

    def expt_pacbio_runs(self, expt):
        """list: PacBio runs for experiment `expt`."""
        assert expt in self.experiments, f"invalid `expt` {expt}"
        return (self.pacbio_df
                .query('experiment == @expt')
                ['pacbio_run']
                .tolist()
                )

    def pacbio_subreads(self, pacbio_run):
        """str: Path to PacBio subreads BAM for `pacbio_run`."""
        assert pacbio_run in self.pacbio_df['pacbio_run'].tolist(), (
                    f"invalid `pacbio_run` {pacbio_run}")
        return str(self.pacbio_df
                   .set_index('pacbio_run')
                   .at[pacbio_run, 'subread_bam']
                   )
