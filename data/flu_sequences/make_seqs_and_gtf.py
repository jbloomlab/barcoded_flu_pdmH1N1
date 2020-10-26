"""Make sequence and GTF files from Genbank plasmid maps."""


import collections
import copy
import glob
import os
import re

import BCBio.GFF

import Bio.Data.IUPACData
import Bio.Seq
import Bio.SeqFeature
import Bio.SeqIO
import Bio.SeqRecord

import yaml


def main():
    """Main body of script."""
    # flu genes
    gene_names = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M', 'NS']

    # viral tag variants
    viral_tags = ['wt', 'syn']
    assert 1 <= len(viral_tags) <= 2, 'script only works for 1 or 2 viral tags'
    tag_file_suffixes = {'wt': '', 'syn': '-dblSyn'}  # suffixes of files

    # map nucleotide sets to their ambiguous characters
    nts_to_ambig_char = {frozenset(vals): char for char, vals in
                         Bio.Data.IUPACData.ambiguous_dna_values.items()}

    # information used to identify features
    polyA = re.compile('A{5}')  # run of at least 5 A nucleotides
    u12 = re.compile('AGC[GA]AAAGCAGG')  # U12 with position 4 polymorphism
    u13 = re.compile('CCTTGTTTCTACT')  # reverse of U13
    # Define introns removed to create NS2 and M2 as in Fig 3 of
    # http://mbio.asm.org/content/5/3/e00070-14.abstract
    mRNA_splice = {
            'M': re.compile('AAC(?P<intron>GTA[CT]GTTC[ACGT]+'
                            'GAAAATTT[AG]CAG)[GA]C'),
            'NS': re.compile('CAG(?P<intron>GTAGA[CT]TG[ACGT]+'
                             'C[CT]TC[TC][TCA]T[TG]CCAG)GA'),
            }
    viralbclen = 16  # length of viral barcode

    genes = []
    tags = collections.defaultdict(dict)

    plasmidmaps = glob.glob('plasmid_maps/*.gb')

    for gene_short_name in gene_names:
        print(f"\nProcessing maps for {gene_short_name}")

        genes_by_tag = {}
        plasmids_by_tag = {}
        for viral_tag in viral_tags:
            suffix = tag_file_suffixes[viral_tag]
            plasmidmatcher = re.compile(f".*_?pH.*{gene_short_name}(_G155E)?"
                                        rf"{suffix}(_\w+)?\.gb")
            plasmidfile = [f for f in plasmidmaps if plasmidmatcher.search(f)]
            if len(plasmidfile) == 1:
                print(f"Reading map for tag {viral_tag} from {plasmidfile[0]}")
                plasmidseq = str(Bio.SeqIO.read(plasmidfile[0], 'genbank').seq)
                plasmids_by_tag[viral_tag] = os.path.basename(plasmidfile[0])
            else:
                raise ValueError(f"not 1 map for {gene_short_name}{viral_tag}")

            # get vRNA flanked by U12 / U13
            u12match = list(u12.finditer(plasmidseq))
            assert len(u12match) == 1, f"Not 1 U12 for {gene_short_name}"
            u13match = list(u13.finditer(plasmidseq))
            assert len(u13match) == 1, f"Not 1 U13 for {gene_short_name}"
            geneseq = plasmidseq[u12match[0].start(0): u13match[0].end(0)]
            print(f"vRNA length is {len(geneseq)}")
            genes_by_tag[viral_tag] = geneseq
        if 1 != len({len(geneseq) for geneseq in genes_by_tag.values()}):
            raise ValueError('all tag variants not the same length')

        geneseq = genes_by_tag[viral_tags[0]]  # annotate this tag version
        plasmid_name = plasmids_by_tag[viral_tags[0]]
        description = (f"influenza {gene_short_name} vRNA from "
                       f"plasmid {plasmid_name}")
        if len(viral_tags) > 1:
            other_plasmids = [plasmids_by_tag[tag] for tag in viral_tags[1:]]
            description += (' with ambiguous nucleotides at sites of viral '
                            f"tags relative to {', '.join(other_plasmids)}")
        gene_name = "flu" + gene_short_name

        gene = Bio.SeqRecord.SeqRecord(
                Bio.Seq.Seq(geneseq),
                id=gene_name,
                description=description,
                name=gene_name,
                annotations={'molecule_type': 'DNA'},
                features=[Bio.SeqFeature.SeqFeature(
                    Bio.SeqFeature.FeatureLocation(0, len(geneseq)),
                    id=gene_name,
                    type='exon',
                    strand=1,
                    qualifiers={
                        "source": "JesseBloom",
                        "gene_id": gene_name,
                        "gene_name": gene_name,
                        "gene_biotype": "vRNA",
                        "transcript_id": gene_name,
                        "label": gene_name,
                        },
                    )],
                )

        # Get mRNA flanked by one nt into u12 and polyA, see:
        # 5' end: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4446424/
        # polyA: https://www.ncbi.nlm.nih.gov/pubmed/7241649
        mrna_start = 1
        polyAmatches = list(polyA.finditer(geneseq))
        assert len(polyAmatches) >= 1, f"no polyA for {gene_short_name}"
        mrna_end = polyAmatches[-1].start(0)
        mrna_name = gene_name + {False: '', True: '1'}[gene_short_name in
                                                       mRNA_splice]
        gene.features.append(
                Bio.SeqFeature.SeqFeature(
                    Bio.SeqFeature.FeatureLocation(mrna_start, mrna_end),
                    id=mrna_name,
                    type='mRNA',
                    strand=1,
                    qualifiers={
                        "source": "JesseBloom",
                        "label": mrna_name,
                        "gene_biotype": "protein_coding",
                        },
                    )
                )

        # is there second splice form?
        if gene_short_name in mRNA_splice:
            intron = list(mRNA_splice[gene_short_name].finditer(geneseq))
            assert len(intron) == 1, f"not 1 intron for {gene_short_name}"
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
                            "source": "JesseBloom",
                            "label": mrna_name,
                            "gene_biotype": "protein_coding",
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
            last_triplet = mrna_start + 3 * (len(mrna_seq[mrna_start:]) // 3)
            protlen = len(mrna_seq[mrna_start: last_triplet].translate(
                    to_stop=True))
            mrna_end = mrna_start + 3 * protlen + 3  # include stop codon
            print("{0} mRNA is of length {1}, encodes protein of {2} residues"
                  .format(mrna.id, len(mrna), protlen))

            # construct CDS location
            assert mrna_start < len(mrna.location.parts[0]), (
                    "start not in first exon")
            assert len(mrna) - mrna_end < len(mrna.location.parts[-1]), (
                    "end not in last exon")
            cds_start = mrna.location.parts[0].start + mrna_start
            cds_end = mrna.location.parts[-1].end - (len(mrna) - mrna_end)
            assert cds_end <= len(gene)
            location_parts = [Bio.SeqFeature.FeatureLocation(
                    cds_start, min(cds_end, mrna.location.parts[0].end))]
            for p in mrna.location.parts[1: -1]:
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
                            "source": "JesseBloom",
                            "label": mrna.id,
                            "gene_biotype": "CDS",
                            "translation": prot
                            },
                        )
                    )

        # annotate viral barcodes in segments that contain them
        if 'N' in geneseq:
            print('Annotating viral barcode...')
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

        # annotate viral tags
        itag = 0
        for i in range(len(geneseq)):
            nts = {tag: seq[i] for tag, seq in genes_by_tag.items()}
            if len(set(nts.values())) == 1:
                continue  # not a tag if all variants the same
            itag += 1
            tagname = f"viral_tag_{itag}"
            tags[gene_name][tagname] = nts
            gene.features.append(
                    Bio.SeqFeature.SeqFeature(
                        Bio.SeqFeature.FeatureLocation(i, i + 1),
                        id=tagname,
                        type=tagname,
                        strand=1,
                        qualifiers=tags[gene_name][tagname],
                        )
                    )
            gene.seq = (
                    gene.seq[:i] +
                    Bio.Seq.Seq(nts_to_ambig_char[frozenset(nts.values())]) +
                    gene.seq[i + 1:]
                    )
        print(f"Annotated {len(tags[gene_name])} tags for {gene_name}")

        genes.append(gene)

    # now write output files
    genefile = 'flu-CA09.fasta'
    print(f"\nWriting {len(genes)} genes to {genefile}")
    Bio.SeqIO.write(genes, genefile, 'fasta')

    gtffile = os.path.splitext(genefile)[0] + '.gtf'
    print(f"Writing gene annotations to {gtffile}")
    just_genes = []
    for gene in genes:
        just_gene = copy.deepcopy(gene)
        just_gene.features = just_gene.features[: 1]
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

    # write Genbank file
    genbankfile = os.path.splitext(genefile)[0] + '.gb'
    print(f"\nWriting annotated genes to {genbankfile}")
    Bio.SeqIO.write(genes, genbankfile, 'genbank')

    # write YAML file with viral tag identities
    tagfile = os.path.splitext(genefile)[0] + '_viral_tags.yaml'
    print(f"\nWriting viral tag identities to {tagfile}")
    with open(tagfile, 'w') as f:
        f.write(yaml.safe_dump(dict(tags)))


# run the script
if __name__ == '__main__':
    main()
