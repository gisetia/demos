# %%
import pandas as pd
from itertools import groupby
from operator import itemgetter

gene = 'CD274'
data_path = '/media/data/nas_scratch/guizela/data'
refseq_dir = f'{data_path}/refseq-processed'


def contig_list_lims(lst):
    '''Divide a list into contiguous sublists and return the lower and
    upper limits of such sublists. (Right open???)
    '''
    lims = []
    for k, g in groupby(enumerate(lst), lambda x: x[0] - x[1]):
        x = list(map(itemgetter(1), g))
        lims.append([x[0], x[-1] + 1])
    return lims


def get_exon_regions(refseq: pd.DataFrame, gene: str) -> pd.DataFrame:
    '''Returns dataframe with start and end positions of all gene exons and
    whether they are cds or utr. Eg.:
    name2	chrom	name	        exon_id	tx_id	region_type	start	end
    CD274	chr9	NM_001267706.1	1	    1	    UTR	        5450502	5450596
    CD274	chr9	NM_001267706.1	2	    1	    UTR	        5456099	5456113
    CD274	chr9	NM_001267706.1	2	    1	    CDS	        5456113	5456165
    '''

    gene_pos = refseq.query('name2 == @gene')

    exon = gene_pos.copy(deep=True) 
    exon = exon.reset_index(drop=True)
    exon['tx_id'] = exon.index + 1

    cols = ['exonStarts', 'exonEnds']
    exon[cols] = exon[cols].applymap(lambda x: x.split(','))

    exon['cdsRange'] = exon.apply(lambda x: range(x.cdsStart, x.cdsEnd),
                                  axis=1)
    exon['exonRange'] = (exon.apply(lambda t:
                                    [range(int(x), int(y)) for i, x
                                     in enumerate(t['exonStarts']) for j, y
                                     in enumerate(t['exonEnds']) if i == j
                                     and x != '' and y != ''], axis=1))

    exon = exon.explode('exonRange')
    exon = exon.reset_index(drop=True)
    exon['exon_id'] = exon.index + 1
    exon = exon.drop(columns=['#bin', 'exonStarts', 'exonEnds', 'score',
                              'cdsStartStat', 'cdsEndStat', 'exonFrames',
                              'cdsStart', 'cdsEnd'])
    exon['exCds_lst'] = exon.apply(lambda t: [x for x in t.exonRange
                                              if x in t.cdsRange], axis=1)
    exon['exUtr_lst'] = exon.apply(lambda t: [x for x in t.exonRange
                                              if x not in t.cdsRange], axis=1)

    exon['CDS'] = exon.apply(lambda t: contig_list_lims(t.exCds_lst), axis=1)
    exon['UTR'] = exon.apply(lambda t: contig_list_lims(t.exUtr_lst), axis=1)

    exon_regions = exon.melt(id_vars=['name2', 'chrom', 'name', 'exon_id',
                                      'tx_id'], value_vars=['CDS', 'UTR'],
                             var_name='region_type', value_name='region_lims')
    exon_regions = exon_regions.explode('region_lims').dropna()
    exon_regions[['start', 'end']] = exon_regions['region_lims'].apply(
        lambda x: pd.Series([x[0], x[1]]))
    exon_regions = exon_regions.drop('region_lims', axis=1)

    exon_regions = exon_regions.sort_values(
        ['tx_id', 'exon_id', 'start']).reset_index(drop=True)

    return exon_regions


refseq = pd.read_parquet(f'{refseq_dir}/ncbi-genes-hg38.parquet.snappy',
                         engine='pyarrow').query('coding == True and '
                                                 'known == True')
exons = get_exon_regions(refseq, gene)
exons

# %%
