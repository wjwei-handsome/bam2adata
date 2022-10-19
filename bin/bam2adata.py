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

## import built-in modules
import argparse
import os
import os.path
import sys

## import local modules
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR + '/src/')

from badmodule import adata
from badmodule.process import batch_process


def get_args():

    # define arguments
    parser = argparse.ArgumentParser(
        description=None, prog='bam2adata.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # required arguments
    required_parser = parser.add_argument_group('required arguments')
    required_parser.add_argument('-g', '--gtf', action='store', dest='in_gtf', type=str, required=True,
                                 help='your gtf file')
    required_parser.add_argument('-b', '--bam', action='store', dest='in_bam', type=str, required=True,
                                 help='your bam file')
    required_parser.add_argument('-d', '--dir', action='store', dest='workdir', type=str, required=True,
                                 help='output dir')
    required_parser.add_argument('--bin', action='store', dest='bin', type=str, required=True,
                                 help='bin size')
    required_parser.add_argument('-p','--pos', action='store', dest='pos_range', type=str, required=True,
                                 help='pos range file path')

    # options arguments
    option_parser = parser.add_argument_group('optional arguments')
    option_parser.add_argument('--indelsize', action='store', dest='indelsize', type=int, default=5,
                               help='allowd indel size to fix the hole')
    option_parser.add_argument('--flanksize', action='store', dest='flanksize', type=int, default=5,
                               help='flanksize')
    option_parser.add_argument('--threads', action='store', dest='threads', type=int, default=4,
                               help='threads')
    option_parser.add_argument('--rmdup', action='store', dest='rmdup', type=bool, default=False,
                               help='remove duplicate reads')
    return parser.parse_args()


if __name__ == '__main__':
    args = get_args()
    ## mkdir  workdir
    try:
        os.mkdir(args.workdir)
        os.mkdir(args.workdir + '/tmp_raw_gem')
    except:
        pass
    # os.mkdir(args.workdir) # TODO check if exists
    # os.mkdir(args.workdir + '/tmp_raw_gem') # TODO make it a tmp dir
    tmp=[int(_) for _ in args.pos_range.split(',')]
    spatial_range = (range(*tmp[0:2]),range(*tmp[2:4]))
    batch_process(args.in_gtf, args.flanksize, args.workdir, args.in_bam, args.indelsize, args.rmdup, args.threads)
    sig = os.system(f"cat {args.workdir}/tmp_raw_gem/*.gem > {args.workdir}/raw_all.gem")
    sig2 = os.system(f"grep -w -f {args.pos_range} {args.workdir}/raw_all.gem > {args.workdir}/raw_gem_pos.gem")
    adata.gem2adata(args.workdir + '/raw_gem_pos.gem', args.bin, args.workdir)
