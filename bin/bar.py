# """

# Demonstrates the use of multiple Progress instances in a single Live display.

# """

# from time import sleep

# from rich.live import Live
# from rich.panel import Panel
# from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn
# from rich.table import Table


# job_progress = Progress(
#     "{task.description}",
#     SpinnerColumn(),
#     BarColumn(),
#     TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
# )
# job1 = job_progress.add_task("[green]Cooking", total=100)
# job2 = job_progress.add_task("[magenta]Baking", total=200)
# job3 = job_progress.add_task("[cyan]Mixing", total=400)

# total = sum(task.total for task in job_progress.tasks)
# overall_progress = Progress()
# overall_task = overall_progress.add_task("All Jobs", total=int(total))

# progress_table = Table.grid()
# progress_table.add_row(
#     Panel.fit(
#         overall_progress, title="Overall Progress", border_style="green", padding=(2, 2)
#     ),
#     Panel.fit(job_progress, title="[b]Jobs", border_style="red", padding=(1, 2)),
# )

# with Live(progress_table, refresh_per_second=20):
#     while not overall_progress.finished:
#         sleep(0.1)
#         for job in job_progress.tasks:
#             if not job.finished:
#                 job_progress.advance(job.id)

#         completed = sum(task.completed for task in job_progress.tasks)
#         overall_progress.update(overall_task, completed=completed)

"""
A rudimentary URL downloader (like wget or curl) to demonstrate Rich progress bars.
"""

import os.path
import sys
from concurrent.futures import ThreadPoolExecutor
import signal
from functools import partial
from threading import Event
from typing import Iterable
from urllib.request import urlopen

from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TaskID,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

progress = Progress(
    TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
    BarColumn(bar_width=None),
    "[progress.percentage]{task.percentage:>3.1f}%",
    "•",
    DownloadColumn(),
    "•",
    TransferSpeedColumn(),
    "•",
    TimeRemainingColumn(),
)


done_event = Event()


def handle_sigint(signum, frame):
    done_event.set()


signal.signal(signal.SIGINT, handle_sigint)


def copy_url(task_id: TaskID, url: str, path: str) -> None:
    """Copy data from a url to a local file."""
    progress.console.log(f"Requesting {url}")
    response = urlopen(url)
    # This will break if the response doesn't contain content length
    progress.update(task_id, total=int(response.info()["Content-length"]))
    with open(path, "wb") as dest_file:
        progress.start_task(task_id)
        for data in iter(partial(response.read, 32768), b""):
            dest_file.write(data)
            progress.update(task_id, advance=len(data))
            if done_event.is_set():
                return
    progress.console.log(f"Downloaded {path}")


def download(urls: Iterable[str], dest_dir: str):
    """Download multiple files to the given directory."""

    with progress:
        with ThreadPoolExecutor(max_workers=4) as pool:
            for url in urls:
                filename = url.split("/")[-1]
                dest_path = os.path.join(dest_dir, filename)
                task_id = progress.add_task("download", filename=filename, start=False)
                pool.submit(copy_url, task_id, url, dest_path)


if __name__ == "__main__":
    # Try with https://releases.ubuntu.com/20.04/ubuntu-20.04.3-desktop-amd64.iso
    if sys.argv[1:]:
        download(sys.argv[1:], "./")
    else:
        print("Usage:\n\tpython downloader.py URL1 URL2 URL3 (etc)")