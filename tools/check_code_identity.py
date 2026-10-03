"""Check that .py files differ only in comments, docstrings and layout.

Compares every .py file under NEW_DIR with the file at the same relative path
under REF_DIR. Both are parsed with ast, which drops comments and layout,
docstrings are removed, and the dumped trees are compared. Use it after an
edit that should change comments, docstrings or spacing only.

With --map, a CSV whose first line is old_path,new_path, each file is first
looked up in the old-to-new name map, so a file that was moved or renamed is
compared with its old version. Paths are relative to REF_DIR and NEW_DIR. An
empty old_path lists a file added on purpose, an empty new_path one removed on
purpose. Rows about other files than .py are ignored, so one table can serve
the MATLAB check too. A map entry naming a .py file that does not exist is an
error.

Every file must be accounted for: a new file with no reference file, and a
reference file that no new file was compared with, are failures, so a moved
file cannot pass uncompared. Two folders with no .py file at all are an error,
since an empty comparison would pass without checking anything. Hidden folders
(git's folder, worktrees, the environments), __pycache__ and venv* are skipped
on both sides.

Prints one line per file, then a summary. Statuses: "same code", and
"added (listed)" or "removed (listed)" for files the map lists; the failures
are "CODE CHANGED", "NO COUNTERPART", "ONLY IN REF" and "SYNTAX ERROR" (a file
that does not parse, or is not UTF-8 text).
Exit code 1 if there is any failure.

    python check_code_identity.py NEW_DIR REF_DIR [--map NAME_MAP.csv]
"""

import argparse
import ast
import csv
import sys
from pathlib import Path

FAILURES = ("CODE CHANGED", "NO COUNTERPART", "ONLY IN REF", "SYNTAX ERROR")
_SKIPPED = ("__pycache__",)


def is_skipped(name):
    """Folders that hold no code of ours: git, worktrees, caches, environments."""
    return name in _SKIPPED or name.startswith((".", "venv"))


def python_files(root):
    """The .py files under `root`, as sorted relative posix paths."""
    if not root.is_dir():
        raise FileNotFoundError(f"folder not found: {root}")
    files = []
    for folder, dirs, names in root.walk():
        # prune in place, so the skipped folders are never walked into
        dirs[:] = [d for d in dirs if not is_skipped(d)]
        for name in names:
            if name.endswith(".py"):
                files.append((folder / name).relative_to(root).as_posix())
    return sorted(files)


def read_name_map(path):
    """The old-to-new name map as a list of (old, new) posix paths."""
    if path is None:
        return []

    # utf-8-sig: a file saved by a spreadsheet can start with a byte order mark
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    if not rows or [cell.strip() for cell in rows[0]] != ["old_path", "new_path"]:
        raise ValueError(f"the first line of {path} must be old_path,new_path")

    pairs = []
    for line_number, row in enumerate(rows[1:], start=2):
        if not any(cell.strip() for cell in row):
            continue
        if len(row) != 2:
            raise ValueError(
                f"line {line_number} of {path} is not old_path,new_path: {row}"
            )
        old, new = (cell.strip().replace("\\", "/") for cell in row)

        # one table can serve every language: keep the rows about .py files
        if old.endswith(".py") or new.endswith(".py"):
            pairs.append((old, new))
    return pairs


def strip_docstrings(tree):
    """Remove the docstring of the module and of every function and class."""
    for node in ast.walk(tree):
        if not isinstance(
            node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ):
            continue

        # a docstring is a string constant as the first statement
        body = node.body
        has_docstring = (
            body
            and isinstance(body[0], ast.Expr)
            and isinstance(body[0].value, ast.Constant)
            and isinstance(body[0].value.value, str)
        )
        if has_docstring:
            # keep the body valid when the docstring was its only statement
            node.body = body[1:] or [ast.Pass()]
    return tree


def code_of(path):
    """Code of a file without comments, docstrings and layout."""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    return ast.dump(strip_docstrings(tree), include_attributes=False)


def compare_pair(new_path, ref_path):
    """Status of one new file against its reference file."""
    # a file that is not UTF-8, or holds a null byte, does not parse either
    try:
        same = code_of(new_path) == code_of(ref_path)
    except (SyntaxError, ValueError):
        return "SYNTAX ERROR"
    if same:
        return "same code"
    return "CODE CHANGED"


def main(new_dir, ref_dir, map_file=None):
    """Compare the two folders, print the result, return the exit code."""
    # the files on both sides, and the name map
    new_dir = Path(new_dir)
    ref_dir = Path(ref_dir)
    new_files = python_files(new_dir)
    ref_files = python_files(ref_dir)
    name_map = read_name_map(map_file)

    # an empty comparison would pass without checking anything
    if not new_files and not ref_files:
        raise FileNotFoundError(f"no .py files to compare in {new_dir} and {ref_dir}")

    # every map entry must name a file that exists: a typo would otherwise leave
    # the file it meant uncompared
    for old, new in name_map:
        if old and old not in ref_files:
            raise FileNotFoundError(
                f"the name map lists {old}, which is not in {ref_dir}"
            )
        if new and new not in new_files:
            raise FileNotFoundError(
                f"the name map lists {new}, which is not in {new_dir}"
            )

    # each new file can come from one old file only
    listed_new = [new for old, new in name_map if new]
    if len(set(listed_new)) < len(listed_new):
        raise ValueError("a new path appears twice in the name map")

    # the old path of each new path in the map, and the files removed on purpose
    old_of = {new: old for old, new in name_map if new}
    listed_removed = {old for old, new in name_map if not new}

    # each new file against its reference file
    rows = []
    compared = set()
    for rel in new_files:
        # the reference file: from the map, or at the same path
        if rel in old_of:
            ref_rel = old_of[rel]
        elif rel in ref_files:
            ref_rel = rel
        else:
            ref_rel = ""

        # added on purpose, no reference file, or compared
        if not ref_rel and rel in old_of:
            status = "added (listed)"
        elif not ref_rel:
            status = "NO COUNTERPART"
        else:
            status = compare_pair(new_dir / rel, ref_dir / ref_rel)
            compared.add(ref_rel)
        rows.append((rel, ref_rel, status))

    # reference files no new file was compared with
    for ref_rel in ref_files:
        if ref_rel in compared:
            continue
        if ref_rel in listed_removed:
            rows.append(("", ref_rel, "removed (listed)"))
        else:
            rows.append(("", ref_rel, "ONLY IN REF"))

    # one line per file, the reference named when it is not at the same path
    counts = {}
    for rel, ref_rel, status in rows:
        counts[status] = counts.get(status, 0) + 1
        if rel == ref_rel:
            print(f"{status:17s} {rel}")
        else:
            print(f"{status:17s} {rel or '-'}  (ref {ref_rel or '-'})")

    # the count of each status; exit code 1 on any failure
    summary = ", ".join(f"{n} {status}" for status, n in counts.items())
    print(f"\n{len(new_files)} new and {len(ref_files)} reference files: {summary}")
    failed = sum(counts.get(status, 0) for status in FAILURES)
    if failed:
        return 1
    return 0


if __name__ == "__main__":
    # the two folders and the name map; the exit code is main's
    parser = argparse.ArgumentParser(description="code identity of two .py folders")
    parser.add_argument("new_dir", help="the code after the edit")
    parser.add_argument("ref_dir", help="the code before the edit")
    parser.add_argument(
        "--map",
        dest="map_file",
        default=None,
        help="CSV old_path,new_path of moved, added and removed files",
    )
    args = parser.parse_args()
    sys.exit(main(args.new_dir, args.ref_dir, args.map_file))
