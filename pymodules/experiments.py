"""Modules for parsing experiments configuration."""


import collections
import copy
import textwrap

import pandas as pd

import ruamel.yaml


class Experiments:
    """Parse and access experiments configuration.

    Parameters
    ----------
    experiments_config : str or dict
        Name of YAML file with experiment configuration, or dict
        read from such a file.

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
    viral_barcodes : dict
        Paths to viral barcode sequencing FASTQ files in dict format
    viral_barcodes_df : pandas.DataFrame
        Data frame with all viral barcode sequencing paths

    """

    def __init__(self, experiments_config):
        """See main class docstring."""

        if isinstance(experiments_config, dict):
            self.config_dict = copy.deepcopy(experiments_config)
        else:
            yaml = ruamel.yaml.YAML(typ='safe')
            with open(experiments_config_file) as f:
                self.config_dict = yaml.load(f)
            assert isinstance(self.config_dict, dict)

        self.experiments = list(self.config_dict)
        self._expect_ncells = {}

        # now get data for each experiment
        valid_keys = {'description',
                      'lab_notes',
                      'expect_ncells',
                      'transcriptomics',
                      'viral_barcodes'}
        transcriptomics_records = []
        viral_barcodes_records = []
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
                        for gene, gene_d in tag_d.items():
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
                                                        ])
        self.transcriptomic_runs = (self.transcriptomics_df
                                    ['transcriptomic_run']
                                    .tolist()
                                    )
        assert (len(self.transcriptomic_runs) ==
                len(set(self.transcriptomic_runs)))
        
        self.viral_barcodes_df = pd.DataFrame(viral_barcodes_records,
                                             columns=['experiment',
                                                      'source',
                                                      'tag',
                                                      'gene',
                                                      'replicate',
                                                      'run',
                                                      'fastq_path'])

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
