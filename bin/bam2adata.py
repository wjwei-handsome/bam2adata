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
import signal
import sys
from concurrent.futures import ThreadPoolExecutor
from itertools import chain
from threading import Event
from typing import Iterable

## import dependencies
import pandas as pd
from rich.progress import (BarColumn, Progress, TaskID, TextColumn,
                           TimeRemainingColumn)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_DIR + '/src/')

from badmodule import adata, gtf
from badmodule.bam import SubBam
from badmodule.utils import get_time


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
                                 help='your gtf file')
    required_parser.add_argument('-d', '--dir', action='store', dest='workdir', type=str, required=True,
                                 help='work dir')
    required_parser.add_argument('--bin', action='store', dest='bin', type=str, required=True,
                                 help='bin size')
    required_parser.add_argument('--pos', action='store', dest='pos_range', type=str, required=True,
                                 help='pos range file path')

    # options arguments
    option_parser = parser.add_argument_group('optional arguments')
    option_parser.add_argument('--indelsize', action='store', dest='indelsize', type=int, default=5,
                               help='allowd indel size to fix the hole')
    option_parser.add_argument('--flanksize', action='store', dest='flanksize', type=int, default=5,
                               help='flanksize')
    option_parser.add_argument('--threads', action='store', dest='threads', type=int, default=1,
                               help='threads')
    option_parser.add_argument('--rmdup', action='store_true', dest='rmdup', default=False,
                               help='remove duplicate reads')
    return parser.parse_args()


progress = Progress(
    TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
    BarColumn(bar_width=None),
    "[progress.percentage]{task.percentage:>3.1f}%",
    "•",
    # DownloadColumn(),
    # "•",
    # TransferSpeedColumn(),
    # "•",
    TimeRemainingColumn(),
)


done_event = Event()


def handle_sigint(signum, frame):
    done_event.set()


signal.signal(signal.SIGINT, handle_sigint)


def _run_bam2rawgem(task_id: TaskID, gene_items:list, path: str) -> None:
    """Copy data from a url to a local file."""
    # progress.console.log(f"Process {task_id} start")
    # response = urlopen(url)
    # This will break if the response doesn't contain content length
    progress.update(task_id, total=len(gene_items))
    # with open(path, "wb") as dest_file:
    progress.start_task(task_id)
    dd_list = []
    for i in gene_items:
        chro, start, end = gtf.get_region(i)
        # gene_id = i['gene_id']
        # start = str(start)
        # end = str(end)
        subregion = SubBam(args.in_bam, chro,start,end,list(chain(*i['ivl'])),args.indelsize)
        d = subregion.bam2rawgem()
        dd_list.append(d)
        progress.update(task_id, advance=1)
        if done_event.is_set():
            return
    merge_df = pd.concat(dd_list)
    merge_df.to_csv(path,sep='\t',index=False)
        # for data in iter(partial(response.read, 32768), b""):
        #     dest_file.write(data)
        #     progress.update(task_id, advance=len(data))
        #     if done_event.is_set():
        #         return
    progress.console.log(f"Done in {path}")


def batch_run(chroms: Iterable[str], dest_dir: str):
    """Download multiple files to the given directory."""

    with progress:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for idx,chrom in  enumerate(chroms):
                # filename = chrom.split("/")[-1]
                # dest_path = os.path.join(dest_dir, filename)
                filename = f"{dest_dir}/batch_{idx}.gem"
                task_id = progress.add_task("download", filename=filename,start=False)
                pool.submit(_run_bam2rawgem, task_id, chrom, filename)

args = get_args()
a = gtf.get_genes(args.in_gtf, args.flanksize)

if __name__ == '__main__':
    batch_size = 5
    batches = [a[i:i + batch_size] for i in range(0, len(a), batch_size)]
    ## mkdir  workdir
    os.mkdir(args.workdir) # TODO check if exists
    os.mkdir(args.workdir + '/tmp_raw_gem') # TODO make it a tmp dir
    batch_run(batches[:5], f'{args.workdir}/tmp_raw_gem')
    sig = os.system(f"cat {args.workdir}/tmp_raw_gem/*.gem > {args.workdir}/raw_all.gem")
    sig2 = os.system(f"grep -f {args.pos_range} {args.workdir}/raw_all.gem > {args.workdir}/raw_gem_pos.gem")
    adata.gem2adata(args.workdir + '/raw_gem_pos.gem', args.bin, args.workdir)


    # @get_time
    # def test(i):
    #     # for i in track(item_list,description='test'):
    #     print('子进程: {} - 任务{}'.format(os.getpid(), i))
    #     chro, start, end = gtf.get_region(i)
    #     gene_id = i['gene_id']
    #     # start = str(start)
    #     # end = str(end)
    #     subregion = SubBam(args.in_bam, chro,start,end,list(chain(*i['ivl'])),args.indelsize)
    #     d = subregion.bam2rawgem(f'data/{gene_id}.outgem')


    # print("CPU内核数:{}".format(cpu_count()))
    # print('当前母进程: {}'.format(os.getpid()))
    # start = time.time()
    # p = Pool(4)
    # for i in a[:10]:
    #     p.apply_async(test, args=(i,))
    # print('等待所有子进程完成。')
    # p.close()
    # p.join()
    # end = time.time()
    # print("总共用时{}秒".format((end - start)))

    # print('当前母进程: {}'.format(os.getpid()))
    # start = time.time()
    # test(a[0])
    # test(a[1])
    # end = time.time()
    # print("总共用时{}秒".format((end - start)))


    # start = time.time()
    # pool = ThreadPoolExecutor(max_workers=8)
    # tasks = []
    # for item in a[:100]:
    #     tasks.append(pool.submit(test, item))
    # for task in track(as_completed(tasks), total=len(tasks)):
    #     if task.result():
    #         raise (task.result())
    # end = time.time()
    # print(f"multi 10 took {(end - start)}")

    # start = time.time()
    # with ThreadPoolExecutor(max_workers=8) as pool:
    # # 使用线程执行map计算
    # # 后面元组有3个元素，因此程序启动3次线程来执行action函数
    #     results = pool.map(test, a[:30])
    # end = time.time()
    # print(f"multi 10 took {(end - start)}")

    # start = time.time()
    # dd_list = []
    # for i in track(a,description='test'):
    #     chro, start, end = gtf.get_region(i)
    #     gene_id = i['gene_id']
    #     # start = str(start)
    #     # end = str(end)
    #     subregion = SubBam(args.in_bam, chro,start,end,list(chain(*i['ivl'])),args.indelsize)
    #     d = subregion.bam2rawgem()
    #     dd_list.append(d)
    # merge_df = pd.concat(dd_list)
    # merge_df.to_csv(args.out_gem,sep='\t',index=False)
    # end = time.time()
    # print(f"no multi 10 took {(end - start)}")
