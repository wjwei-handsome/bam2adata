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
from rich.progress import track
import pandas as pd


HIT_DICT = {0: 'EXONIC', 1: 'INTRONIC', 2: 'INTERGENIC'}


class Bam:

    # gene_tag = "GE:Z"
    # hit_tag = "XF:i"
    # spatial_x_tag = "Cx:i"
    # spatial_y_tag = "Cy:i"

    def __init__(self, bamfile, exon_ivl, indel_size,):
        self.bamfile = pysam.AlignmentFile(bamfile, 'rb')
        self.exon_ivl = exon_ivl
        self.indel_size = indel_size
        self.contig_reads = self.get_total_reads(self.bamfile)
        # self.multi_blocks = 0
        # self.splice_beyond = 0
        # self.intergeneic = 0
        # self.unmapgene = 0
        # self.total_reads = self.get_total_reads(bamfile)



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




class BamContig:
    gene_tag = "GE:Z"
    hit_tag = "XF:i"
    spatial_x_tag = "Cx:i"
    spatial_y_tag = "Cy:i"

    def __init__(self, bamfile, contig, exon_ivl, indel_size, total_reads):
        self.bamfile = bamfile
        self.exon_ivl = exon_ivl
        self.indel_size = indel_size
        self.multi_blocks = 0
        self.splice_beyond = 0
        self.intergeneic = 0
        self.unmapgene = 0
        self.contig_name = contig
        self.bam_contig = self.bamfile.fetch(contig=contig)
        self.total_reads = total_reads
        # self.rawgem_df = pd.DataFrame(columns=['gene','spatial','MIDCount','status','chrtype','segments'])
        # self.total_reads = self.get_total_reads(bamfile)


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

    def bam2rawgem(self, outrawgem: str):
        with open(outrawgem, 'w') as out_gem_handler:
            for read in track(
                    self.bam_contig,
                    description=f"[bold green]Procesing bam {self.contig_name}",
                    total=self.total_reads
                ):
                print(self.total_reads)
                status = 'UN'
                if read.is_unmapped or read.is_qcfail:
                    continue
                else:
                    refname = read.reference_name
                    chrtype = 'Chrom' if refname.startswith('chr') else 'Contig'
                    try:
                        gene = read.get_tag(self.gene_tag)
                        hit = read.get_tag(self.hit_tag)
                        x_pos = read.get_tag(self.spatial_x_tag)
                        y_pos = read.get_tag(self.spatial_y_tag)
                        spatial = f'{x_pos}_{y_pos}'
                        pos = read.reference_start + 1
                        segments = self.parse_cigar(
                            read.cigartuples, pos, read.is_reverse)
                        if HIT_DICT[hit] == 'INTONIC':
                            status = 'Unspliced'
                        elif HIT_DICT[hit] == 'EXONIC':
                            if 'N' in read.cigarstring:
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

                        writeline = '\t'.join(
                            [gene, spatial, '1', status, chrtype, str(segments)])

                        out_gem_handler.write(writeline + '\n')

                    except KeyError:
                        self.unmapgene += 1
                        continue
            self.stat_df = self.get_stats()

    def get_stats(self):
        stat_df = pd.DataFrame(
            [{
            'total_reads': self.total_reads,
            'unmapgene': self.unmapgene,
            'intergeneic': self.intergeneic,
            'splice_beyond': self.splice_beyond,
            'multi_blocks': self.multi_blocks
            }]
            )
        return stat_df

pd.DataFrame()