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

from badmodule import gtf
from badmodule.bam import SubBam

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


def _run_bam2rawgem(task_id: TaskID, gene_items:list, path: str, bamfile_path: str, indelsize: int, rmdup: bool) -> None:
    """run bam2rawgem for a batch of genes"""
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
        subregion = SubBam(bamfile_path, chro,start,end,list(chain(*i['ivl'])),indelsize)
        d = subregion.bam2rawgem(rmdup)
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


def batch_run(chroms: Iterable[str], dest_dir: str, bamfile: str, indelsize: int, rmdup: bool, threads: int) -> None:
    """Download multiple files to the given directory."""

    with progress:
        with ThreadPoolExecutor(max_workers=threads) as pool:
            for idx,chrom in  enumerate(chroms):
                # filename = chrom.split("/")[-1]
                # dest_path = os.path.join(dest_dir, filename)
                filename = f"{dest_dir}/batch_{idx}.gem"
                task_id = progress.add_task("batch_run", filename=filename,start=False)
                pool.submit(_run_bam2rawgem, task_id, chrom, filename, bamfile, indelsize, rmdup)


def batch_process(in_gtf, flank_size, workdir, bamfile, indelsize, rmdup, threads):
    a = gtf.get_genes(in_gtf, flank_size)

    ## simulate 5 batches
    batch_size = 1000
    batches = [a[i:i + batch_size] for i in range(0, len(a), batch_size)]
    batch_run(batches, f'{workdir}/tmp_raw_gem', bamfile, indelsize, rmdup, threads)
