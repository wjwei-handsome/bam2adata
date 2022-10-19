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

from itertools import chain

import numpy as np
import pandas as pd
import pysam
from rich.progress import track

HIT_DICT = {0: 'EXONIC', 1: 'INTRONIC', 2: 'INTERGENIC'}

POS_FILE = 'data/pos.range'
pos_range = pd.read_csv(POS_FILE, header=None, names=[
                                'pos', ]).values.flatten()

class Bam:

    # gene_tag = "GE:Z"
    # hit_tag = "XF:i"
    # spatial_x_tag = "Cx:i"
    # spatial_y_tag = "Cy:i"

    def __init__(self, bamfile, exon_ivl, indel_size,):
        self.bamfile = pysam.AlignmentFile(bamfile, 'rb')
        self.indel_size = indel_size
        self.contig_reads = self.get_total_reads(self.bamfile)




    @staticmethod
    def get_total_reads(bamfile):
        try:
            idx_stats = bamfile.get_index_statistics()
            contig_reads = {x.contig: x.total for x in idx_stats}
            # TODO add logger
        except ValueError:
            pass # TODO: add logger;check idx file?

        return contig_reads

    def split_contig_bam(self, scaf: bool = False):

        contigs = self.bamfile.references
        if not scaf:
            contigs = [x for x in contigs if not x.startswith('scaf')]

        contig_bams = [
            BamContig(
                self.bamfile,
                x,
                self.exon_ivl,
                self.indel_size,
                self.contig_reads[x]
                ) for x in contigs
            ]
        return contig_bams




class SubBam:
    gene_tag = "GE:Z"
    hit_tag = "XF:i"
    spatial_x_tag = "Cx:i"
    spatial_y_tag = "Cy:i"

    def __init__(self, bamfile, chro,start,end, exon_ivl, indel_size):
        self.bamfile = pysam.AlignmentFile(bamfile, 'rb')
        self.exon_ivl = exon_ivl
        self.indel_size = indel_size
        self.multi_blocks = 0
        self.splice_beyond = 0
        self.intergeneic = 0
        self.unmapgene = 0
        self.region = f"{chro}:{start}-{end}"
        self.bam_region = self.bamfile.fetch(contig=chro,start=start,end=end)
        # self.total_reads = total_reads
        # self.rawgem_df = pd.DataFrame(columns=['gene','spatial','MIDCount','status','chrtype','segments'])
        # self.total_reads = self.get_total_reads(bamfile)

    @staticmethod
    def pos_bin_mapper(x, binWidth=50):
        """trans pos('18183_12175') to bin('18150_12150'), default binWidth 50

        Args:
            x ([str in pandas.colnums]): [Any]
            binWidth ([int]): [the width of bin, default 50]

        Returns:
            [str]: [x and y concat with '_' after bin]
        """
        x_num_str_before_bin = x.split('_')[0]
        y_num_str_before_bin = x.split('_')[1]
        x_num_int_after_bin = np.floor(
            int(x_num_str_before_bin)/binWidth).astype(int)*binWidth
        y_num_int_after_bin = np.floor(
            int(y_num_str_before_bin)/binWidth).astype(int)*binWidth
        new_pos = str(x_num_int_after_bin)+'_'+str(y_num_int_after_bin)
        return new_pos


    def parse_cigar(self, cigartuples: tuple, pos: int, reverse: bool) -> list:
        """_summary_

        Args:
            cigar (tuple): cigar string
            pos (int): start position
            reverse (bool): if reverse strand
        Returns:
            list: segments list
        """
        segments = []
        hole_to_remove = set()
        ref_skip = False
        clip5 = clip3 = 0
        p = pos
        cigartuples.reverse() if reverse else None
        for i, (operation_id, length) in enumerate(cigartuples):
            if operation_id == 0:  # vcy.CIGAR[operation_id] == "BAM_CMATCH"
                segments.append((p, p + length - 1))
                p += length
            # A splice || vcy.CIGAR[operation_id] == 'BAM_CREF_SKIP'
            elif operation_id == 3:
                ref_skip = True
                p += length
            # A deletion || cy.CIGAR[operation_id] == 'BAM_CDEL'
            elif operation_id == 2:
                if length <= self.indel_size:
                    try:
                        if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                            hole_to_remove.add(len(segments) - 1)
                    except IndexError:
                        pass
                p += length
            # bases at 5' or 3' are NOT part of the alignment || vcy.CIGAR[operation_id] == 'BAM_CSOFT_CLIP'
            elif operation_id == 4:
                if p == pos:
                    clip5 = length  # At start of alignment
                else:
                    # Must be at end of alignment vcy.CIGAR[operation_id] in ["BAM_CINS", "BAM_CHARD_CLIP"]
                    clip3 = length
                p += length
            elif operation_id == 1:  # An insertion BAM_CINS
                if length <= self.indel_size:
                    try:
                        if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                            hole_to_remove.add(len(segments) - 1)
                    except IndexError:
                        pass
                # else do nothing
                # NOTE: maybe we should make so that the reads get discarded
            elif operation_id == 5:  # BAM_CHARD_CLIP
                print("Hard clip was encountered! All mapping are assumed soft clipped")

        # Merge segments separated by small insertions and deletions
        # NOTE maybe sorted is not required realy
        for a, b in enumerate(sorted(hole_to_remove)):
            segments[b - a] = (segments.pop(b - a)[0], segments[b - a][1])

        return segments

    def bam2rawgem(self,rmdup: bool):
        aaa = []
        for read in self.bam_region:
            status = 'UN'
            if rmdup:
                if read.is_duplicate:
                    continue
            if read.is_unmapped or read.is_qcfail:
                self.unmapgene += 1
                continue
            else:
                # refname = read.reference_name
                # chrtype = 'Chrom' if refname.startswith('chr') else 'Contig'
                try:
                    gene = read.get_tag(self.gene_tag)
                    hit = read.get_tag(self.hit_tag)
                    x_pos = read.get_tag(self.spatial_x_tag)
                    y_pos = read.get_tag(self.spatial_y_tag)
                    spatial = f'{x_pos}_{y_pos}'
                    # if spatial not in pos_range:
                    #     continue
                    pos = read.reference_start + 1
                    if HIT_DICT[hit] == 'INTRONIC':
                        status = 'Unspliced'
                    elif HIT_DICT[hit] == 'EXONIC':
                        if 'N' in read.cigarstring:
                            segments = self.parse_cigar(read.cigartuples, pos, read.is_reverse)
                            starts, ends = zip(*segments)
                            # if segemnts number >=4, then it is multi-blocks
                            if len(segments) >= 4:
                                self.multi_blocks += 1
                                continue
                            elif len(segments) == 3:
                                if any([_ in self.exon_ivl for _ in (starts[1:3] + ends[:2])]):
                                    # if starts[1] in self.exon_ivl and starts[2] in self.exon_ivl and ends[0] in self.exon_ivl and ends[1] in self.exon_ivl:
                                    status = 'Spliced'
                                else:
                                    self.splice_beyond += 1
                                    continue
                            else:
                                if any([_ in self.exon_ivl for _ in [starts[1], ends[0]]]):
                                    # if starts[1] in self.exon_ivl and ends[0] in self.exon_ivl:
                                    status = 'Spliced'
                                else:
                                    self.splice_beyond += 1
                                    continue
                        else:
                            status = 'Ambiguous'
                    elif HIT_DICT[hit] == 'INTERGENIC':
                        self.intergeneic += 1
                        continue
                    else:
                        continue
                except KeyError:
                    self.intergeneic += 1
                    continue
            aaa.append([gene, spatial, 1, status])
        aaa_df = pd.DataFrame(aaa, columns=['gene', 'spatial', 'MIDCount', 'status'])
        aaa_df['spatial'] = aaa_df['spatial'].map(self.pos_bin_mapper)
        aaaa = aaa_df.groupby(['gene', 'spatial', 'status']).size().reset_index().rename(columns={0: 'MIDcounts'})
        # self.stat_df = self.get_stats()
        return aaaa


    def get_stats(self):
        stat_df = pd.DataFrame(
            [{
            # 'total_reads': self.total_reads,
            'unmapgene': self.unmapgene,
            'intergeneic': self.intergeneic,
            'splice_beyond': self.splice_beyond,
            'multi_blocks': self.multi_blocks
            }]
            )
        return stat_df
