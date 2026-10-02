"""Compare two output folders file by file.

For a change that must not change results: run the old and the new code on
separate copies of the inputs, then compare what they wrote. Every file is
compared by type:

    identical bytes        same
    .csv .tsv              same table (pandas), else DIFFERENT
    .npy .npz              same arrays, dtypes and shapes, NaN equal to NaN
    .png .jpg              same pixels; "same render" if at most 1 grey level
                           differs on a few pixels (renderer noise)
    .tif .tiff             same pixels on every page of a stack, exactly: these
                           hold data, where one grey level is a real difference
    .eps                   same text without the creation date and title lines
    .txt .json .md .py .m  same text
    anything else          DIFFERENT if the bytes differ

Also lists the files present in only one folder; a file that cannot be read is
"compare failed", and the others are still compared. MATLAB's .mat and .fig
files and videos are compared by bytes here, and they hold the time they were
written: compare the MATLAB pipelines' folders, and videos, with
tools/sep_compare_outputs.m, which reads their contents.

--replace OLD NEW replaces text in the old run's tables and text files before
comparing, for the paths a run writes, which name its own tree. --newer-than
reports a file of NEW_DIR last written before that time as "NOT REWRITTEN", so
a run that wrote nothing cannot pass; compare output folders only with it,
since inputs are never rewritten.

Needs numpy, pandas and pillow. Exit code 1 if any file is not the same.

    python compare_outputs.py OLD_DIR NEW_DIR [--ignore REGEX ...]
        [--replace OLD NEW] [--newer-than "2026-10-01 09:00"]
"""

import argparse
import filecmp
import io
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

SAME = ("same", "same render")


def read_text(path, replacements):
    """Text of a file, with the (old, new) replacements applied.

    Read as latin-1, which maps every byte to one character, so no byte is
    lost or merged with another whatever the file's encoding.
    """
    text = path.read_text(encoding="latin-1")
    for old, new in replacements:
        text = text.replace(old, new)
    return text


def same_table(a, b, replacements):
    """Same table, whatever the files' line endings and number formatting."""
    sep = ","
    if a.suffix.lower() == ".tsv":
        sep = "\t"
    table_a = pd.read_csv(io.StringIO(read_text(a, replacements)), sep=sep)
    table_b = pd.read_csv(io.StringIO(read_text(b, [])), sep=sep)
    if table_a.equals(table_b):
        return "same"
    return "DIFFERENT"


def load_arrays(path):
    """The arrays of a .npy or .npz file, by name."""
    if path.suffix.lower() == ".npy":
        return {"x": np.load(path, allow_pickle=False)}
    with np.load(path, allow_pickle=False) as npz:
        return {name: npz[name] for name in npz.files}


def same_arrays(a, b):
    """Same arrays, dtypes and shapes, NaN equal to NaN."""
    arrays_a = load_arrays(a)
    arrays_b = load_arrays(b)
    if sorted(arrays_a) != sorted(arrays_b):
        return "DIFFERENT"

    # NaN counts as equal to NaN only where the dtype can hold it
    for name, x in arrays_a.items():
        y = arrays_b[name]
        if x.dtype != y.dtype or x.shape != y.shape:
            return "DIFFERENT"
        if not np.array_equal(x, y, equal_nan=x.dtype.kind in "fc"):
            return "DIFFERENT"
    return "same"


def same_image(a, b, render_tolerance):
    """Same pixels on every frame of two images.

    With `render_tolerance`, "same render" when at most 1 grey level differs,
    on under 0.1% of the pixels of a frame.
    """
    status = "same"
    with Image.open(a) as image_a, Image.open(b) as image_b:
        n_frames = getattr(image_a, "n_frames", 1)
        if getattr(image_b, "n_frames", 1) != n_frames:
            return "DIFFERENT"

        # one frame at a time, so a large stack is never held twice in memory
        for index in range(n_frames):
            image_a.seek(index)
            image_b.seek(index)
            pixels_a = np.asarray(image_a).astype(np.int64)
            pixels_b = np.asarray(image_b).astype(np.int64)
            if pixels_a.shape != pixels_b.shape:
                return "DIFFERENT"
            diff = np.abs(pixels_a - pixels_b)
            if diff.size == 0 or diff.max() == 0:
                continue

            # renderer noise: tiny, on a small fraction of the pixels
            if render_tolerance and diff.max() <= 1 and (diff > 0).mean() < 0.001:
                status = "same render"
                continue
            return "DIFFERENT"
    return status


def strip_eps_dates(text):
    """EPS text without the lines that change at every export."""
    return re.sub(r"%%(CreationDate|Title):[^\n]*", "", text)


def compare(a, b, replacements):
    """Status of one pair of files."""
    if filecmp.cmp(a, b, shallow=False):
        return "same"

    # bytes that differ: compare by type, as the module docstring lists
    suffix = a.suffix.lower()
    if suffix in (".csv", ".tsv"):
        return same_table(a, b, replacements)
    if suffix in (".npy", ".npz"):
        return same_arrays(a, b)
    if suffix in (".png", ".jpg", ".jpeg"):
        return same_image(a, b, render_tolerance=True)
    if suffix in (".tif", ".tiff"):
        return same_image(a, b, render_tolerance=False)
    if suffix == ".eps":
        text_a = strip_eps_dates(read_text(a, replacements))
        text_b = strip_eps_dates(read_text(b, []))
    elif suffix in (".txt", ".json", ".md", ".py", ".m"):
        text_a = read_text(a, replacements)
        text_b = read_text(b, [])
    else:
        return "DIFFERENT"
    if text_a == text_b:
        return "same"
    return "DIFFERENT"


def relative_files(root, skip):
    """Relative posix paths of the files under `root`, without the skipped ones."""
    if not root.is_dir():
        raise FileNotFoundError(f"folder not found: {root}")
    files = set()
    for path in root.rglob("*"):
        rel = path.relative_to(root).as_posix()
        if path.is_file() and not any(pattern.search(rel) for pattern in skip):
            files.add(rel)
    return files


def main(old_dir, new_dir, ignore, replacements, newer_than=None):
    """Compare the two folders, print the result, return the exit code."""
    # the files on both sides, without the ignored ones
    old_dir = Path(old_dir)
    new_dir = Path(new_dir)
    skip = [re.compile(pattern) for pattern in ignore]
    old_files = relative_files(old_dir, skip)
    new_files = relative_files(new_dir, skip)

    # an empty comparison would pass without checking anything
    if not old_files and not new_files:
        raise FileNotFoundError(f"no files to compare in {old_dir} and {new_dir}")

    # each file of either folder, in one folder only or compared by type; a line
    # for each file that is not the same
    counts = {}
    for rel in sorted(old_files | new_files):
        detail = ""
        if rel not in new_files:
            status = "only in old"
        elif rel not in old_files:
            status = "only in new"
        else:
            try:
                status = compare(old_dir / rel, new_dir / rel, replacements)
            except (OSError, ValueError) as err:
                # a file that cannot be read is reported, and the others still compared
                status = "compare failed"
                detail = f"  ({err})"

            # a file the run did not write can only repeat what was there before
            if newer_than is not None:
                written = datetime.fromtimestamp((new_dir / rel).stat().st_mtime)
                if written < newer_than:
                    detail = (
                        f"  (written {written:%Y-%m-%d %H:%M:%S}, compared: {status})"
                    )
                    status = "NOT REWRITTEN"

        counts[status] = counts.get(status, 0) + 1
        if status != "same":
            print(f"{status:14s} {rel}{detail}")

    # the count of each status; exit code 1 if any file is not the same
    print("\n" + ", ".join(f"{n} {status}" for status, n in sorted(counts.items())))
    n_bad = sum(n for status, n in counts.items() if status not in SAME)
    if n_bad:
        return 1
    return 0


if __name__ == "__main__":
    # the two folders and the options; the exit code is main's
    parser = argparse.ArgumentParser(description="compare two output folders")
    parser.add_argument("old_dir", help="outputs of the reference run")
    parser.add_argument("new_dir", help="outputs of the run after the change")
    parser.add_argument(
        "--ignore",
        nargs="*",
        default=[r"\.log$"],
        help="regular expressions of relative paths to skip",
    )
    parser.add_argument(
        "--replace",
        nargs=2,
        action="append",
        default=[],
        metavar=("OLD", "NEW"),
        help="text replaced in the old run's tables and text files",
    )
    parser.add_argument(
        "--newer-than",
        type=datetime.fromisoformat,
        default=None,
        help="files of NEW_DIR written before this are NOT REWRITTEN",
    )
    args = parser.parse_args()
    sys.exit(main(args.old_dir, args.new_dir, args.ignore, args.replace, args.newer_than))
