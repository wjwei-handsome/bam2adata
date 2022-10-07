import pandas as pd
from tqdm import tqdm
from itertools import chain
class ExonIvl:

    def __init__(self, indel_size: int) -> None:
        self.indel_size = indel_size

def get_geneid_attr(attr: str) -> str:
    return attr.split(';')[0].split(' ')[1].replace('"','')

def get_ivl(x, indel_size):
    return range(x['start']-indel_size,x['end']+indel_size+1)

def get_genes(gtf_file: str, indel_size: int) -> list:

    tqdm.pandas(desc='pandas bar')
    gtf_df = pd.read_csv(gtf_file, sep='\t', header=None, comment='#',usecols=[0,2,3,4,8])
    gtf_df.columns=['chr','type','start','end','attr']
    gtf_df['gene_id'] = gtf_df['attr'].progress_apply(get_geneid_attr)
    del gtf_df['attr']
    exon_df = gtf_df[gtf_df['type'] == 'exon']
    del exon_df['type']
    del gtf_df
    exon_df['ivl'] = exon_df.progress_apply(get_ivl, axis=1, indel_size=indel_size)
    gene_exon_df = exon_df.groupby('gene_id')['ivl'].apply(list).reset_index()
    gene_exon_ivls_dict = pd.Series(gene_exon_df['ivl'].values, index=gene_exon_df['gene_id']).to_dict()
    gene_exon_ivls_dict = {k: chain(*v) for k, v in gene_exon_ivls_dict.items()}

    return gene_exon_ivls_dict

# df['gene_id']=df['attr'].map(lambda x:x.split(';')[0].split(' ')[1].strip('"'))
# del df['attr']
# exon_df = df[df['type'] == 'exon']