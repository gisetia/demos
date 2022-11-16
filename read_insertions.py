# %%
import pandas as pd
from typing import Optional

screen_type = 'ip'  # 'ip', 'pa'
screen_name = 'PDL1_IFNg'

select_by = 'gene'  # Either 'gene' or 'position'
gene = 'CD274'  # Gene symbol in refseq, eg 'CD274'
position = 'chr9:5,450,503-5,470,567' # Position, eg 'chr9:5,450,503-5,470,567'

data_path = '/media/data/nas_scratch/guizela/data'

# Data directories
if screen_type == 'ip':
    ins_dir = f'{data_path}/screen-insertions'
elif screen_type == 'pa':
    ins_dir = f'{data_path}/screen-insertions_activating'
refseq_dir = f'{data_path}/refseq-processed'
#%%

def read_insertions(data_dir: str, screen_name: str, chrom: str, start: int,
                    end: int, assembly: str = 'hg38', trim_length: int = 50,
                    screen_type: str = 'ip') -> pd.DataFrame:

    data_path = (f'{data_dir}/{screen_name}/'
                 f'{assembly}/{trim_length}/')

    insertions = None

    if screen_type == 'ip' or screen_type == 'pa':
        filename = 'insertions.parquet.snappy'
        try:
            insertions = pd.read_parquet(f'{data_path}{filename}')
            insertions = insertions.query('chr == @chrom '
                                          '& pos >= @start '
                                          '& pos <= @end')
        except FileNotFoundError:
            print(f'!! Error: no file found for screen {screen_name} of '
                  f'type {screen_type}')

    return insertions


def read_refseq(data_dir: str, assembly: str = 'hg38',
                coding_only: Optional[bool] = True,
                known_only: Optional[bool] = True) -> pd.DataFrame:

    filename = f'{data_dir}/ncbi-genes-{assembly}.parquet.snappy'
    refseq = pd.read_parquet(filename, engine='pyarrow')

    if coding_only:
        refseq = refseq[refseq['coding']]
    if known_only:
        refseq = refseq[refseq['known']]

    return refseq


def convert_position(position: str) -> tuple:

    chrom = f'chr{position.split(":")[0][3:]}'

    # Convert to 0-based left-closed right-open
    start = int(position.split(':')[1].split('-')[0].replace(',', '')) - 1
    end = int(position.split(':')[1].split('-')[1].replace(',', ''))

    return (chrom, start, end)


def get_gene_position(gene: str, refseq: pd.DataFrame) -> tuple:

    gene_pos = refseq.query('name2 == @gene')

    chrom = gene_pos.chrom.head(1).values[0]
    start = gene_pos.txStart.min()
    end = gene_pos.txEnd.max()

    strand = gene_pos.strand.head(1).values[0]

    return (chrom, start, end, strand)


def get_insertions(data_dir: str, screen_name: str, select_by: str = 'gene',
                   gene: str = None, position: str = None) -> pd.DataFrame:

    if select_by == 'gene':
        if gene is None:
            raise KeyboardInterrupt('!! Error: no gene name was provided')
        else:
            print(f'Getting insertions for gene {gene} in screen '
                  f'{screen_name}')
    elif select_by == 'position':
        if position is None:
            raise KeyboardInterrupt('!! Error: no position was provided')
        else:
            print(f'Getting insertions for interval {position} in screen '
                  f'{screen_name}')

    chrom = start = end = None
    if select_by == 'gene':
        refseq = read_refseq(refseq_dir)
        chrom, start, end, strand = get_gene_position(gene, refseq)
        print(f'Gene {gene} is on strand {strand}')
    elif select_by == 'position':
        chrom, start, end = convert_position(position)

    insertions = read_insertions(data_dir, screen_name, chrom, start, end,
                                 screen_type=screen_type)

    return insertions


insertions = get_insertions(ins_dir, screen_name, select_by=select_by,
                            gene=gene, position=position)

insertions

# %%
