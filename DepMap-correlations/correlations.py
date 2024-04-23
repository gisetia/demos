# %%
import numpy as np
import pandas as pd
import multiprocessing as mp
from itertools import combinations
from functools import partial
from scipy.stats import pearsonr
from itertools import combinations
from statsmodels.stats.multitest import fdrcorrection

cores = 30


def correlation(i, dataset):

    if i[0] % 1000 == 0 and i[1] % 1000 == 0:
        print(i)

    # ignore nans
    nan = np.logical_or(np.isnan(dataset[i[0]]), np.isnan(dataset[i[1]]))
    corr = pearsonr(dataset[i[0]][~nan], dataset[i[1]][~nan])

    return corr + (*i, len(dataset[i[0]][~nan]))


# %%
# CRISPR gene effect correlations

# crispr = pd.read_parquet('gs://gisetia-ccle/processed_data/cell_cols/'
#                          'DepMap23Q4/CRISPRGeneEffect.pq')

# gene_pairs = combinations(range(len(crispr)), 2)
# partial_func = partial(correlation, dataset=crispr.values)


# with mp.Pool(cores) as p:
#     result = p.imap(partial_func, gene_pairs, chunksize=100000)
#     result = list(result)

# crispr_corr = pd.DataFrame(result, columns=['corr', 'p', 'gene1',
#                                             'gene2', 'n'])
# crispr_corr['p_fdr'] = fdrcorrection(crispr_corr['p'])[1]

# crispr_corr[['gene1', 'gene2']] = crispr_corr[['gene1', 'gene2']].applymap(
#     lambda x: crispr.index[x])

# crispr_corr.to_parquet('processed_data/crispr_corr.pq')

# %%
# Expression correlations

exp = pd.read_parquet('gs://gisetia-ccle/processed_data/cell_cols/DepMap23Q4'
                      '/OmicsExpressionProteinCodingGenesTPMLogp1.pq')

gene_pairs = combinations(range(len(exp)), 2)
partial_func = partial(correlation, dataset=exp.values)

with mp.Pool(cores) as p:
    result = p.imap_unordered(partial_func, gene_pairs, chunksize=100000)
    result = list(result)

exp_corr = pd.DataFrame(result, columns=['corr', 'p', 'gene1',
                                         'gene2', 'n'])
exp_corr['p_fdr'] = fdrcorrection(exp_corr['p'])[1]

exp_corr[['gene1', 'gene2']] = exp_corr[['gene1', 'gene2']].applymap(
    lambda x: exp.index[x])
exp_corr.to_parquet('processed_data/exp_corr.pq')

# # %%
# # Merge data

# corr = pd.merge(crispr_corr, exp_corr, left_on=['gene1', 'gene2'],
#                 right_on=['gene1', 'gene2'], suffixes=('_crispr', '_exp'),
#                 how='outer').query('p_fdr_crispr < 0.05 or p_fdr_exp < 0.05')
# corr.to_parquet('processed_data/corr.pq')


# %%
