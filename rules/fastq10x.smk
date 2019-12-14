"""``snakemake`` rules related to generating 10X FASTQ files."""


rule fastq10x_qc_stats:
    """Aggregates and plots the QC stats for the 10X FASTQ files."""
    input:
        qc_stats=expand(join(config['fastq10x_dir'], "{run10x}_qc_stats.csv"),
                        run10x=illumina_runs_10x.index)
    output:
        qc_stats=report(join(config['fastq10x_dir'], 'fastq10x_qc_stats.csv'),
                        caption='../report/fastq10x_qc_stats.rst',
                        category='10X FASTQ files',
                        ),
        qc_plot=report(join(config['fastq10x_dir'], 'fastq10x_qc_stats.svg'),
                       caption='../report/fastq10x_qc_plot.rst',
                       category='10X FASTQ files',
                       )
    run:
        # read and aggregate the stats
        print('Reading 10X FASTQ QC stats from:\n\t' +
              '\n\t'.join(input.qc_stats))
        stats = pd.concat([(pd.read_csv(statfile,
                                        names=['statistic', 'value'])
                            .assign(run10x=run10x)
                            )
                           for statfile, run10x in zip(input.qc_stats,
                                                       illumina_runs_10x.index)
                           ],
                          ignore_index=True)
        print(f"Writing the aggregated stats to {output.qc_stats}")
        stats.to_csv(output.qc_stats, index=False)

        # plot the stats
        _ = (ggplot(stats, aes('run10x', 'value')) +
             geom_point(size=2) +
             facet_wrap('~ statistic', ncol=4, scales='free_y') +
             theme(axis_text_x=element_text(angle=90),
                   figure_size=(12, 4), panel_spacing_x=0.6) +
                   expand_limits(y=(0, 1)) +
             scale_y_continuous(
                labels=mizani.formatters.custom_format('{:.2g}'),
                )
             ).save(output.qc_plot, verbose=False, bbox_inches='tight')
                           

rule make_fastq10x:
    """Make 10X FASTQ files from BCL runs using `cellranger mkfastq`."""
    output:
        fastqR1=join(config['fastq10x_dir'], "{run10x}_all_R1.fastq.gz"),
        fastqR2=join(config['fastq10x_dir'], "{run10x}_all_R2.fastq.gz"),
        mkfastq10x_dir=directory(join(config['mkfastq10x_dir'], "{run10x}")),
        qc_stats=join(config['fastq10x_dir'], "{run10x}_qc_stats.csv"),
        csv=temp("_mkfastq_{run10x}.csv"),  # input for cellranger mkfastq
        mro=temp("__{run10x}.mro"),  # created by cellranger mkfastq
    params:
        run10x="{run10x}"
    threads:
        config['max_cpus']
    run:
        # write CSV file for `cellranger mkfastq`
        with open(output.csv, 'w') as f:
            f.write('Lane,Sample,Index\n' +
                    ','.join(map(str,
                             [illumina_runs_10x.at[params.run10x, 'lane'],
                              params.run10x,
                              illumina_runs_10x.at[params.run10x, 'index']])))

        # run `cellranger mkfastq`
        cmds = ['cellranger', 'mkfastq',
                '--run', illumina_runs_10x.at[params.run10x, 'bcl_folder'],
                '--id', params.run10x,  # output directory name
                '--csv', output.csv,
                '--delete-undetermined',
                '--qc',
                f"--localcores={threads}",
                ]
        print(f"\nRunning the following commands:\n{' '.join(cmds)}\n")
        subprocess.check_call(cmds)

        # move `cellranger mkfastq` output to desired location
        print(f"\nMoving `cellranger mkfastq` output from {params.run10x} "
              f"to {output.mkfastq10x_dir}\n")
        shutil.move(params.run10x, output.mkfastq10x_dir)

        # get names of R1 and R2 FASTQ files from `cellranger mkfastq` output
        fastq_glob = f"{output.mkfastq10x_dir}/outs/fastq_path/*/*/*.fastq.gz"
        fastqs = sorted(glob.glob(fastq_glob))
        fastqregex = re.compile(r'_(?P<read>R1|R2|I1)_\d{3}\.fastq\.gz')
        r1s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R1']
        r2s = [f for f in fastqs if fastqregex.search(f).group('read') == 'R2']
        print(f"\nFASTQ files matching {fastq_glob}:\n"
              f"  R1 files: {','.join(map(basename, r1s))}\n"
              f"  R2 files: {','.join(map(basename, r2s))}\n")
        assert len(r1s) == len(r2s)

        # concatenate all R1 and R2 files into merged FASTQs for run
        print(f"\nCreating merged FASTQs {output.fastqR1}, {output.fastqR2}\n")
        for outfq, fqlist in [(output.fastqR1, r1s), (output.fastqR2, r2s)]:
            with open(outfq, 'wb') as f:
                subprocess.call(['cat'] + fqlist, stdout=f)

        # extract sample QC stats from `cellranger mkfastq` JSON into CSV
        qc_json = join(config['mkfastq10x_dir'], params.run10x,
                       'outs/qc_summary.json')
        print(f"\nExtracting QC stats from {qc_json} to {output.qc_stats}\n")
        with open(qc_json) as f:
            (pd.Series(json.load(f)['sample_qc'][params.run10x]['all'])
             .to_csv(output.qc_stats, header=False)
             )
