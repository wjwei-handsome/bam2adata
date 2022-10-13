#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :gem2adata.py
@说明        : trans gem to anndata
@时间        :2022/03/18 23:15:06
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''

import argparse
import logging
import math

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse as sp

from badmodule.utils import _error, _info, _warn


def gemZoomer(gem_file):
    """zoom the gem to max(full) size

    Args:
        gem_file (str): gem file path
        clu_list (str): cluster list

    Returns:
        gem: zommed gem
    """

    try:
        # read_type_dict = {"gene": str, "x": int,
        #                   "y": int, "MIDcounts": int, "status": str}
        raw_gem = pd.read_csv(gem_file, sep='\t', header=None)
        raw_gem.columns = ['gene','spatial','status','MIDcounts']
        raw_gem['x'] = raw_gem['spatial'].map(lambda x: int(x.split('_')[0]))
        raw_gem['y'] = raw_gem['spatial'].map(lambda x: int(x.split('_')[1]))
        del raw_gem['spatial']
        gem = raw_gem

    except pd.errors.ParserError:
        _error('Reading gem file Error, Please check !')

    if 'status' not in gem.columns:
        gem['status'] = 'S'
    min_x = gem['x'].min()
    max_x = gem['x'].max()
    min_y = gem['y'].min()
    max_y = gem['y'].max()
    x_scal = max_x-min_x+1
    y_scal = max_y-min_y+1
    _info(f"max x:{max_x} max y: {max_y} ")

    region = (0, x_scal-1, 0, y_scal-1, 'ALL')

    x_scal = region[1]-region[0]+1
    y_scal = region[3]-region[2]+1
    gem['x'] = gem['x']-min_x-region[0]
    gem['y'] = gem['y']-min_y-region[2]
    gem = gem[(gem.x >= 0) & (gem.x < x_scal) &
              (gem.y >= 0) & (gem.y < y_scal)]

    return gem

def transGemToAnnData(gem, bin):
    """trans gem data to adata

    Args:
        gem (gem): gem
        bin (int): bin size
        sample (str): sample name

    Returns:
        adata: adata
    """

    _info("Trans gem to AnnData format.")
    # count info
    bin = int(bin)
    x_pixel = math.ceil((gem['x'].max()-gem['x'].min()+1)/bin)
    y_pixel = math.ceil((gem['y'].max()-gem['y'].min()+1)/bin)
    _info(f'The compressed matrix shape is x={x_pixel},y={y_pixel}')

    # compress the matrix by dataframe
    gem['cx'] = np.floor(gem['x']/bin).astype(int)
    gem['cy'] = np.floor(gem['y']/bin).astype(int)
    gem['bin'] = gem['cx'].astype(str) + '_' + gem['cy'].astype(str)

    # extract the bin and gene categorical in dataframe
    bin_cat = pd.Categorical(gem['bin'])
    '''
    May be bugs here !!!
    '''
    # bin_list = np.sort(np.array(bin_cat.unique()))
    bin_list = np.sort(gem['bin'].unique())  # WWJ
    gem['bin_index'] = bin_cat.codes

    gene_cat = pd.Categorical(gem['gene'])
    # gene_list = np.sort(np.array(gene_cat.unique()))
    gene_list = np.sort(gem['gene'].unique())  # WWJ
    gem['gene_index'] = gene_cat.codes

    # use coo-matrix to compress the gem quickly
    mat = sp.coo_matrix(
        (gem['MIDcounts'], (gem['bin_index'], gem['gene_index']))).tocsr()  # 构建一个稀疏矩阵

    # extract the unspliced mRNA matrix
    un_gem = gem[gem['status'] == 'Unspliced']
    un_counts = np.append(un_gem['MIDcounts'], 0)
    un_bin = np.append(un_gem['bin_index'], gem['bin_index'].max())
    un_gene = np.append(un_gem['gene_index'], gem['gene_index'].max())
    un_mat = sp.coo_matrix((un_counts, (un_bin, un_gene))).tocsr()
    del un_gem

    # extract the spliced mRNA matrix
    sp_gem = gem[gem['status'].isin(['Spliced', 'Ambiguous'])]
    sp_counts = np.append(sp_gem['MIDcounts'], 0)
    sp_bin = np.append(sp_gem['bin_index'], gem['bin_index'].max())
    sp_gene = np.append(sp_gem['gene_index'], gem['gene_index'].max())
    sp_mat = sp.coo_matrix((sp_counts, (sp_bin, sp_gene))).tocsr()
    del sp_gem

    # set the obs list, which contain the sample infos
    obs = pd.DataFrame(index=bin_list)
    obs['bin_loc'] = bin_list
    loc = obs['bin_loc'].str.split('_', 1, expand=True).astype(int)
    loc.columns = ['cx', 'cy']
    obs = pd.merge(obs, loc, how='left', left_index=True, right_index=True)
    obs['sample'] = 'sample'
    # set the var list, which contain the gene info
    var = pd.DataFrame(index=gene_list)
    var['gene_ID'] = gene_list
    # trans the compressed data to AnnData and save as h5ad file
    adata = ad.AnnData(mat, obs=obs, var=var, layers={
                       'spliced': sp_mat, 'unspliced': un_mat})
    # adata = ad.AnnData(mat,obs=obs,var=var,) ##wwj
    return adata



def gem2adata(in_gem, bin, out_dir):
    """main workflow

    Args:
        args (args): args

    Returns:
        adata: adata
        sample: sample name
    """

    # sample = in_gem.split('/')[-1].replace('.gem', '')
    # sample = 'wwjtest'
    gem = gemZoomer(in_gem)
    adata = transGemToAnnData(gem, bin)

    sc.pp.filter_cells(adata, min_counts=1)
    adata.uns['gem_type'] = 'bin'
    adata.obsm['X_spatial'] = np.concatenate((adata.obs.cx, adata.obs.cy), axis=0).\
        reshape(2, -1).T.astype(np.float32)


    out_data = f'{out_dir}/bin_{bin}.h5ad'
    adata.write(out_data, as_dense='X')
