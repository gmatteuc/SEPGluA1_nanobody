"""The numbers of the ISH line for the text: one table per step, then all of them.

Every step of the ISH line writes the numbers it computed, the ones the text quotes,
to tables/numbers_<step>.csv: one row per number, its name, its value as the step
rounded it and what it is. run_ish_overview gathers them into one table with the
step each comes from, so docs/ISH_ANALYSIS.md quotes every number from one file:

    numbers_for_the_text.csv    step, name, value, what
    numbers_for_the_text.txt    the same as aligned lines, for reading

Imported by the modules that write a step's numbers, and by ish.overview.
"""

from pathlib import Path

import pandas as pd

from sepmap.structures import TABLES

ALL_NUMBERS = TABLES / "numbers_for_the_text.csv"
ALL_NUMBERS_TXT = TABLES / "numbers_for_the_text.txt"

# the columns of a step's numbers
COLUMNS = ["name", "value", "what"]


def numbers_path(step: str) -> Path:
    """Where a step writes its numbers: tables/numbers_<step>.csv."""
    return TABLES / f"numbers_{step}.csv"


def numbers_frame(rows: list[tuple]) -> pd.DataFrame:
    """A step's numbers as a table from (name, value, what) rows, values kept as given."""
    return pd.DataFrame(rows, columns=COLUMNS, dtype=object)


def gather_numbers(tables: Path | None = None) -> pd.DataFrame:
    """Every numbers_<step>.csv of the run in one table: step, name, value, what.

    The values stay as the steps wrote them (text), so nothing is rounded twice.
    """
    if tables is None:
        tables = TABLES
    parts = []
    for path in sorted(tables.glob("numbers_*.csv")):
        if path == ALL_NUMBERS:
            continue
        part = pd.read_csv(path, dtype=str, keep_default_na=False)
        part.insert(0, "step", path.stem.removeprefix("numbers_"))
        parts.append(part)
    if not parts:
        raise FileNotFoundError(f"no numbers_*.csv in {tables}: run steps 13 to 29 first")
    return pd.concat(parts, ignore_index=True)


def lookup(numbers: pd.DataFrame) -> dict[str, str]:
    """{'step.name': value}; a name written twice by one step keeps its first value."""
    out = {}
    for step, name, value in zip(numbers["step"], numbers["name"], numbers["value"]):
        out.setdefault(f"{step}.{name}", value)
    return out


def text_lines(numbers: pd.DataFrame) -> list[str]:
    """numbers_for_the_text.txt: one aligned line per number, step by step."""
    width = max(len(f"{s}.{n}") for s, n in zip(numbers["step"], numbers["name"]))
    lines = []
    for step, part in numbers.groupby("step", sort=False):
        lines.append(f"[{step}]")
        for _, r in part.iterrows():
            key = f"{step}.{r['name']}"
            lines.append(f"  {key:<{width}}  {r['value']:<12}  {r['what']}")
        lines.append("")
    return lines


def ordinal(value: str | float) -> str:
    """A rank as words write it: '1st', '2nd', '11th', '17th'."""
    k = int(float(value))
    suffix = "th"
    if k % 100 not in (11, 12, 13):
        suffix = {1: "st", 2: "nd", 3: "rd"}.get(k % 10, "th")
    return f"{k}{suffix}"
