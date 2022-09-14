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

    # options arguments
    option_parser = parser.add_argument_group('optional arguments')
    option_parser.add_argument('--indelsize', action='store', dest='indelsize', type=int, default=5,
                               help='allowd indel size to fix the hole')
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    print(args)
