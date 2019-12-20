"""Make sequence and GTF files from Genbank plasmid maps."""


import re
import os
import glob
import itertools
import collections
import copy

import Bio.Alphabet
import Bio.Data.IUPACData
import Bio.Seq
import Bio.SeqFeature
import Bio.SeqIO
import Bio.SeqRecord

import BCBio.GFF


def main():
    """Main body of script."""

    # flu genes
    gene_names = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

    # plus / minus synonymous viral tags
    viral_tags = ['', '-dblSyn']
    tag_names = {'': 'wt', '-dblSyn': 'syn'}

    # map nucleotide sets to their ambiguous characters
    nts_to_ambig_char = {frozenset(vals): char for char, vals in
                         Bio.Data.IUPACData.ambiguous_dna_values.items()}

    # information used to identify features
    polyA = re.compile('A{5}') # run of at least 5 A nucleotides
    u12 = re.compile('AGC[GA]AAAGCAGG') # U12 with position 4 polymorphism
    u13 = re.compile('CCTTGTTTCTACT') # reverse of U13
    # Define introns removed to create NS2 and M2 as in Fig 3 of
    # http://mbio.asm.org/content/5/3/e00070-14.abstract
    mRNA_splice = {
            'M':re.compile('AAC(?P<intron>GTA[CT]GTTC[ACGT]+GAAAATTT[AG]CAG)[GA]C'),
            'NS':re.compile('CAG(?P<intron>GTAGA[CT]TG[ACGT]+C[CT]TC[TC][TCA]T[TG]CCAG)GA'),
            }
    viralbclen = 16  # length of viral barcode

    genes = collections.defaultdict(list)

    plasmidmaps = glob.glob('plasmid_maps/*.gb')

    for syntype, gene_short_name in itertools.product(viral_tags, gene_names):

        print("\nProcessing {0}{1}...".format(gene_short_name, syntype))

        # read plasmid
        plasmidmatcher = re.compile(f".*_?pH.*{gene_short_name}(_G155E)?"
                                    f"{syntype}(_\w+)?\.gb")
        plasmidfile = [f for f in plasmidmaps if plasmidmatcher.search(f)]

        assert len(plasmidfile) == 1, \
                f"{len(plasmidfile)} maps for {gene_short_name}{syntype}"
        plasmidfile = plasmidfile[0]
        print(f"Reading sequence from {plasmidfile}")
        plasmid = str(Bio.SeqIO.read(plasmidfile, 'genbank').seq)

        # get vRNA flanked by U12 / U13
        u12matches = list(u12.finditer(plasmid))
        assert len(u12matches) == 1, "Not exactly one U12 for {0}{1}".format(
                gene_short_name, syntype)
        u13matches = list(u13.finditer(plasmid))
        assert len(u13matches) == 1, "Not exactly one U13 for {0}{1}".format(
                gene_short_name, syntype)
        geneseq = plasmid[u12matches[0].start(0) : u13matches[0].end(0)]
        print("Gene length is {0}".format(len(geneseq)))

        gene_name = "flu" + gene_short_name

        gene = Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(geneseq, Bio.Alphabet.generic_dna),
                id=gene_name,
                name=gene_name,
                description='{0} vRNA from plasmid {1}'.format(
                        gene_short_name, plasmidfile),
                features=[Bio.SeqFeature.SeqFeature(
                    Bio.SeqFeature.FeatureLocation(0, len(geneseq)),
                    id=gene_name,
                    type='exon',
                    strand=1,
                    qualifiers={
                        "source": "JesseBloom",
                        "gene_id":gene_name,
                        "gene_name":gene_name,
                        "gene_biotype":"vRNA",
                        "transcript_id":gene_name,
                        "label":gene_name,
                        },
                    )],
                )

        # Get mRNA flanked by one nt into u12 and polyA, see:
        # 5' end described: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4446424/
        # polyA described: https://www.ncbi.nlm.nih.gov/pubmed/7241649
        mrna_start = 1
        polyAmatches = list(polyA.finditer(geneseq))
        assert len(polyAmatches) >= 1, "Not more than one polyA for {0}{1}".format(
                gene_short_name, syntype)
        mrna_end = polyAmatches[-1].start(0)
        mrna_name = gene_name + {False:"", True:"1"}[gene_short_name in mRNA_splice]
        gene.features.append(
                Bio.SeqFeature.SeqFeature(
                    Bio.SeqFeature.FeatureLocation(mrna_start, mrna_end),
                    id=mrna_name,
                    type='mRNA',
                    strand=1,
                    qualifiers={
                        "source":"JesseBloom",
                        "label":mrna_name,
                        "gene_biotype":"protein_coding",
                        },
                    )
                )

        # is there second splice form?
        if gene_short_name in mRNA_splice:
            intron = list(mRNA_splice[gene_short_name].finditer(geneseq))
            assert len(intron) == 1, "not exactly one intron for {0}{1}".format(
                    gene_short_name, syntype)
            intron_start = intron[0].start('intron')
            intron_end = intron[0].end('intron')
            mrna_name = gene_name + '2'
            gene.features.append(
                    Bio.SeqFeature.SeqFeature(
                        Bio.SeqFeature.CompoundLocation([
                            Bio.SeqFeature.FeatureLocation(
                                mrna_start, intron_start),
                            Bio.SeqFeature.FeatureLocation(
                                intron_end, mrna_end)
                            ]),
                        id=mrna_name,
                        type="mRNA",
                        strand=1,
                        qualifiers={
                            "source":"JesseBloom",
                            "label":mrna_name,
                            "gene_biotype":"protein_coding",
                            },
                        )
                    )

        # for each mRNA, get the CDS, which we heuristically identify as
        # first ATG initiated sequence
        for mrna in filter(lambda gene: gene.type == 'mRNA', gene.features):

            # find start / end of first ATG initiated sequence
            mrna_seq = mrna.extract(gene).seq
            assert len(mrna_seq) == len(mrna)
            mrna_start = str(mrna_seq).index('ATG')
            if re.match('fluNP', mrna.id):
                # get main, not upstream, ORF for NP:
                # https://www.biorxiv.org/content/biorxiv/early/2019/09/04/738427
                assert mrna_seq[mrna_start + 18: mrna_start + 21] == 'ATG'
                mrna_start = mrna_start + 18
            last_triplet = mrna_start + 3 * (len(mrna_seq[mrna_start : ]) // 3)
            protlen = len(mrna_seq[mrna_start : last_triplet].translate(
                    to_stop=True))
            mrna_end = mrna_start + 3 * protlen + 3 # include stop codon
            print("{0} mRNA is of length {1}, encodes protein of {2} residues"
                    .format(mrna.id, len(mrna), protlen))

            # construct CDS location
            assert mrna_start < len(mrna.location.parts[0]),\
                    "start not in first exon"
            assert len(mrna) - mrna_end < len(mrna.location.parts[-1]),\
                    "end not in last exon"
            cds_start = mrna.location.parts[0].start + mrna_start
            cds_end = mrna.location.parts[-1].end - (len(mrna) - mrna_end)
            assert cds_end <= len(gene)
            location_parts = [Bio.SeqFeature.FeatureLocation(
                    cds_start, min(cds_end, mrna.location.parts[0].end))]
            for p in mrna.location.parts[1 : -1]:
                location_parts.append(p)
            if len(mrna.location.parts) > 1:
                location_parts.append(Bio.SeqFeature.FeatureLocation(
                        mrna.location.parts[-1].start, cds_end))

            if len(location_parts) > 1:
                cds_location = Bio.SeqFeature.CompoundLocation(location_parts)
            else:
                cds_location = location_parts[0]
            prot = cds_location.extract(gene).seq.translate()
            assert len(prot) == protlen + 1
            assert prot[0] == 'M' and prot[-1] == '*'
            gene.features.append(
                    Bio.SeqFeature.SeqFeature(
                        cds_location,
                        id=mrna.id,
                        type="CDS",
                        strand=1,
                        qualifiers={
                            "source":"JesseBloom",
                            "label":mrna.id,
                            "gene_biotype":"CDS",
                            "translation":prot
                            },
                        )
                    )

        # annotate viral barcodes in segments that contain them
        if 'N' in geneseq:
            print(f"Annotating viral barcode...")
            m = re.fullmatch(f"[ACGT]+(?P<viralbc>N{{{viralbclen}}})[ACGT]+",
                             geneseq)
            if not m:
                raise ValueError('could not match viral barcode')
            gene.features.append(
                    Bio.SeqFeature.SeqFeature(
                        Bio.SeqFeature.FeatureLocation(m.start('viralbc'),
                                                       m.end('viralbc')),
                        id='viral_barcode',
                        type='viral_barcode',
                        strand=1
                        )
                    )

        genes[syntype].append(gene)

    # now write FASTA / GTF for wt viral tag variant
    for syntype in viral_tags[:1]:

        genefile = 'flu-CA09{0}.fasta'.format(syntype)
        print(f"\nWriting {len(genes[syntype])} genes to {genefile}")
        Bio.SeqIO.write(genes[syntype], genefile, 'fasta')

        mrnas = [Bio.SeqRecord.SeqRecord(mrna.extract(gene).seq,
                    id=mrna.id, description='')
                for gene in genes[syntype]
                for mrna in gene.features if mrna.type == 'mRNA']
        mrnafile = 'flu-CA09{0}-mRNA.fasta'.format(syntype)
        print("Writing the {0} mRNAs to {1}".format(len(mrnas), mrnafile))
        Bio.SeqIO.write(mrnas, mrnafile, 'fasta')

        gtffile = os.path.splitext(genefile)[0] + '.gtf'
        print("Writing gene annotations to {0}".format(gtffile))
        just_genes = []
        for gene in genes[syntype]:
            just_gene = copy.deepcopy(gene)
            just_gene.features = just_gene.features[ : 1]
            just_genes.append(just_gene)
        with open(gtffile, 'w') as f:
            BCBio.GFF.write(just_genes, f)
        # Hacky conversion of GFF3 file to GTF, changing 9th column to
        # delimit qualifiers with spaces / quotes rather than equals sign.
        # Conversion is probably not robust to all qualifiers, but works here.
        with open(gtffile) as f:
            lines = f.readlines()
        newlines = []
        for line in lines:
            if line[0] == '#':
                newlines.append(line)
            else:
                entries = line.strip().split('\t')
                assert len(entries) == 9, str(len(entries)) + '\n' + line
                newqualifiers = []
                for qualifier in entries[-1].split(';'):
                    (key, value) = qualifier.split('=')
                    newqualifiers.append('{0} "{1}"'.format(key, value))
                entries[-1] = '; '.join(newqualifiers)
                newlines.append('\t'.join(entries) + '\n')
        with open(gtffile, 'w') as f:
            f.write(''.join(newlines))

    # write Genbank annotating the tags
    genbankfile = 'flu-CA09.gb'
    print(f"\nWriting tag-annotated genes to {genbankfile}")
    if len(viral_tags) != 2:
        raise RuntimeError('code for tagged Genbank only works for 2 tags')
    genbank_genes = []
    tags = collections.defaultdict(dict)
    for igene, gene_name in enumerate(gene_names):
        gene = copy.deepcopy(genes[viral_tags[0]][igene])
        gene2 = genes[viral_tags[1]][igene]
        assert gene.features[0].type == 'exon'
        gene.features[0].type = 'vRNA_rc'  # type from exon -> vRNA_rc
        geneseq = str(gene.seq)
        gene2seq = str(gene2.seq)
        assert len(geneseq) == len(gene2seq)
        itag = 0
        for i, (nt, nt2) in enumerate(zip(geneseq, gene2seq)):
            if nt == nt2:
                continue  # not a tag if both variants the same
            itag += 1
            tagname = f"tag_{itag}"
            tags[gene_name][tagname] = {tag_names[viral_tags[0]]: nt,
                                        tag_names[viral_tags[1]]: nt2}
            gene.features.append(
                    Bio.SeqFeature.SeqFeature(
                        Bio.SeqFeature.FeatureLocation(i, i + 1),
                        id=tagname,
                        type=tagname,
                        strand=1,
                        qualifiers=tags[gene_name][tagname],
                        )
                    )
            gene.seq = (gene.seq[:i] +
                        Bio.Seq.Seq(nts_to_ambig_char[frozenset([nt, nt2])],
                                    alphabet=gene.seq.alphabet) +
                        gene.seq[i + 1:])
        print(f"Annotated {len(tags[gene_name])} nucleotide tags for {gene_name}")
        genbank_genes.append(gene)
    Bio.SeqIO.write(genbank_genes, genbankfile, 'genbank')


# run the script
if __name__ == '__main__':
    main()

