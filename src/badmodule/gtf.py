from itertools import chain

import pandas as pd
from tqdm import tqdm

from badmodule.logger import console
from badmodule.utils import _info, _warn


class ExonIvl:

    def __init__(self, flank_size: int) -> None:
        self.flank_size = flank_size

def get_geneid_attr(attr: str) -> str:
    return attr.split(';')[0].split(' ')[1].replace('"','')

def get_ivl(x, flank_size):
    return range(x['start']-flank_size,x['end']+flank_size+1)

def get_genes(gtf_file: str, flank_size: int) -> list:

    with console.status("[bold green]Reading GTF file...") as status:
        # _info('READING...')
        tqdm.pandas()
        gtf_df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#',usecols=[0,2,3,4,8])
        gtf_df.columns=['chr','type','start','end','attr']
        _info("[bold green]READ GTF Done!")
    with console.status("[bold green]Processing GTF file...") as status:
        gtf_df['gene_id'] = gtf_df['attr'].apply(get_geneid_attr)
        del gtf_df['attr']
        exon_df = gtf_df[gtf_df['type'] == 'exon']
        del exon_df['type']
        del gtf_df
        exon_df['ivl'] = exon_df.apply(get_ivl, axis=1, flank_size=flank_size)
        del exon_df['start']
        del exon_df['end']
        gene_exon_df = exon_df.groupby(['gene_id','chr'])['ivl'].apply(list).reset_index()
        # gene_exon_ivls_dict = pd.Series(gene_exon_df['ivl'].values, index=gene_exon_df['gene_id']).to_dict()
        # gene_exon_ivls_dict = {k: chain(*v) for k, v in gene_exon_ivls_dict.items()}
        gene_exon_ivls_dict = gene_exon_df.to_dict('records')
        _info("[bold green]PROCESS GTF Done!")

    return gene_exon_ivls_dict


def get_region(gene_exon_idvls_dict: dict,) -> tuple:
    """_summary_

    Args:
        gene_exon_idvls_dict (dict): gene exon ivls dict

    Returns:
        str: region
    """
    chro = "chr"+str(gene_exon_idvls_dict['chr'])
    # chain_ivls = chain(*gene_exon_idvls_dict['ivl']) # NOTETHIS will cause memory error
    start = min(chain(*gene_exon_idvls_dict['ivl']))
    end = max(chain(*gene_exon_idvls_dict['ivl']))
    return (chro, start, end)
