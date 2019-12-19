"""``snakemake`` rules related to generating 10X FASTQ files."""


rule fastq10x_qc_analysis:
    """Quality control analysis of the 10X FASTQ files."""
    input:
        qc_stats=expand(join(config['fastq10x_dir'], "{run10x}_qc_stats.csv"),
                        run10x=illumina_runs_10x.index),
        nb='notebooks/fastq10x_qc_analysis.ipynb'
    output:
        nb=join(config['fastq10x_dir'], 'fastq10x_qc_analysis.ipynb'),
        nb_html=report(join(config['fastq10x_dir'],
                            'fastq10x_qc_analysis.html'),
                       caption='../report/fastq10x_qc_analysis.rst',
                       category='10X FASTQ files')
    run:
        run_nb_to_html(
                input_nb=input.nb,
                output_nb=output.nb,
                parameters={
                    'illumina_runs_10x': illumina_runs_10x.index.tolist(),
                    'input_qc_stats': input.qc_stats,
                    },
                )


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
