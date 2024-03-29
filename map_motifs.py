import pandas as pd
import subprocess as sp

motifs = ['CTTTAAA', 'CTTTAAG']
save_dir = '.'

bowtie_idx = ('/media/data/reference/genomes/BOWTIE_GRCh38/'
              'GCA_000001405.15_GRCh38_no_alt_analysis_set')
for motif in motifs:
    cmd = (f'bowtie {bowtie_idx} -a -v 0 -c {motif} > '
           f'{save_dir}/{motif}_bowtie.tsv')
    sp.run(cmd, shell=True)

chroms = [str(x) for x in range(1, 23)] + ['X', 'Y']
bed_list = list()
for motif in motifs:
    bowtie = pd.read_csv(f'{save_dir}/{motif}_bowtie.tsv', sep='\t',
                         header=None).rename(columns={2: 'chrom', 1: 'strand',
                                                      4: 'motif'})

    bed = bowtie[['chrom', 'strand', 'motif']]
    bed['chrom'] = bed['chrom'].str.replace('chr', '')
    bed['start'] = bowtie[3]
    bed['end'] = bowtie[3] + len(motif)
    bed['name'] = bowtie.apply(
        lambda x: f'{x["motif"]}-{x["chrom"]}:{x[3]+1}-{x[3]+len(motif)}',
        axis=1)
    bed['score'] = 1

    bed = bed[['chrom', 'start', 'end', 'name', 'score', 'strand']].query(
        'chrom in @chroms')

    bed.sort_values(['chrom', 'start']).to_csv(f'{save_dir}/{motif}.bed',
                                               sep='\t', index=False, 
                                               header=False)
