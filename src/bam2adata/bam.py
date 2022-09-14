#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :bam.py
@说明        :
@时间        :2022/09/14 16:35:09
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''

import pysam


def process_bamfile(bamfile, exon_ivl, indel_size, outgemfile):
    """process bamfile and output gem file

    Args:
        bamfile (str): bamfile path
        exon_ivl (list): flank added exon interval
        indel_size (int): allowed indel size to fix hole
        outgemfile (str): gemfile path

    Returns:
        4 statsis
    """
    multi_blocks = 0
    splice_beyond = 0
    intergeneic = 0
    unmapgene = 0
    bamhandler = pysam.AlignmentFile(bamfile, 'rb')
    with open(outgemfile, 'w') as out_gem_handler:
        for read in tqdm(bamhandler, desc='Process bamfile:'):
            status = 'UN'
            if read.is_unmapped or read.flag >= 256:
                continue
            else:
                refName = read.reference_name
                chrType = 'NORMAL' if refName.startswith('chr') else 'Contig'
                try:
                    gene = read.get_tag("GE:Z")
                    hit = read.get_tag("XF:Z")
                    spatial_pos = read.get_tag("CB:Z")
                    pos = read.reference_start + 1
                    segments = parse_cigar_tuple(
                        read.cigartuples, pos, reverse=read.is_reverse, INDEL_SIZE=indel_size)
                    if hit == 'INTRONIC':
                        status = 'Unspliced'
                    elif hit == 'EXONIC':
                        if 'N' in read.cigarstring:
                            starts, ends = zip(*segments)
                            if len(segments) >= 4:
                                multi_blocks += 1
                                continue
                            elif len(segments) == 3:
                                if starts[1] in exon_ivl and starts[2] in exon_ivl and ends[0] in exon_ivl and ends[1] in exon_ivl:
                                    status = 'Spliced'
                                else:
                                    splice_beyond += 1
                                    continue
                            else:
                                if starts[1] in exon_ivl and ends[0] in exon_ivl:
                                    status = 'Spliced'
                                else:
                                    splice_beyond += 1
                                    continue
                        else:
                            status = 'Ambiguous'
                    else:
                        intergeneic += 1
                        continue

                    writeline = '\t'.join(
                        [gene, spatial_pos, '1', status, chrType, str(segments)])
                    out_gem_handler.writelines(writeline + '\n')
                except:
                    unmapgene += 1
                    continue
    return multi_blocks, splice_beyond, intergeneic, unmapgene
