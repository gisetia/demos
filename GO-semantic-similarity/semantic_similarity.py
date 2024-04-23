#%%


def semantic_similarity(pair, ont, gene_terms, count=[0]):
    # as defined in https://www.nature.com/articles/s41467-019-13058-9#Sec9

    count[0] += 1
    if not count[0] % 50000:
        print(count)

    if gene_terms.get(pair[0]) and gene_terms.get(pair[1]):
        common_terms = (set(gene_terms[pair[0]])
                        .intersection(gene_terms[pair[1]]))
    else:
        common_terms = {}

    if common_terms == {}:
        sem_sim = 0
        min_term = 'None'
        min_num = np.inf
        min_name = 'None'
    else:
        common_ont = (ont.query('term_id in @common_terms')
                      .sort_values('size').iloc[0])

        min_term = common_ont['term_id']
        min_num = common_ont['size']
        sem_sim = 2/min_num
        min_name = common_ont['name']

    return *pair, sem_sim, min_term, min_num, min_name
