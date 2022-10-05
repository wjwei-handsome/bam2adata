#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :bam2adata.py
@说明        :
@时间        :2022/09/14 16:15:13
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''

import argparse
from rich.progress import track

import sys
from typing import *
import os
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR + '/src/')

from badmodule.bam import Bam


def get_args():
    # define arguments
    parser = argparse.ArgumentParser(
        description=None, prog='bam2adata.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('--gtf', action='store', dest='in_gtf', type=str, required=True,
                                 help='your gtf file')
    required_parser.add_argument('--bam', action='store', dest='in_bam', type=str, required=True,
                                 help='your gtf file')
    required_parser.add_argument('--gem', action='store', dest='out_gem', type=str, required=True,
                                 help='output gem file name')
    # required_parser.add_argument('--dir', action='store', dest='dir', type=str, required=True,
    #                              help='tmp work dir')

    # options arguments
    option_parser = parser.add_argument_group('optional arguments')
    option_parser.add_argument('--indelsize', action='store', dest='indelsize', type=int, default=5,
                               help='allowd indel size to fix the hole')
    return parser.parse_args()

def add_pos_flank(ivl, pos_str, size):
    """add flank to a interval list

    Args:
        ivl (list): interval list
        pos_str (str): postion
        size (int): flank size
    """
    pos = int(pos_str)
    for i in range(pos-size, pos+size+1):
        ivl.append(i)


def get_exon_ivl(gtf_file):
    """get exon interval list

    Args:
        gtf_file (str): gtf file path

    Returns:
        list: exon interval list
    """
    exon_ivl = []
    with open(gtf_file, 'r') as gtf_file_handler:
        for line in track(gtf_file_handler, description='Process gtf file'):
            if not line.startswith('#'):
                fields = line.rstrip().split('\t')
                chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
                if 'exon' in feature_type:
                    add_pos_flank(exon_ivl, start_str, 5)
                    add_pos_flank(exon_ivl, end_str, 5)


    return exon_ivl


if __name__ == '__main__':
    args = get_args()
    exon_ivl = get_exon_ivl(args.in_gtf)
    bam = Bam(args.in_bam, exon_ivl, args.indelsize)
    bamcontigs = bam.split_contig_bam()
    test_chr10 = bamcontigs[-1]
    test_chr10.bam2rawgem(args.out_gem)
    # print(args)
