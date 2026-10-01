"""Where the project's data lives, and the settings of the Python route.

The data root is worked out the way get_paths.m works it out for MATLAB: the
parent of the code folder, plus data. The code folder is the one that holds
get_paths.m, found by walking up from this file, so the answer does not
depend on how deep in the code folder this file sits (mapping/sepmap/).

SEP_DATA_ROOT points the whole data tree somewhere else. The refactor's checks
run on copies of the data and must never write into the real one, and two
guards, the same as in get_paths.m, make that hard to get wrong: a copy of the
code may not use the production data, and the snapshot on G: is a backup,
never a data root.

settings.toml, next to the package in mapping/, is read once into SETTINGS.
print_settings gives every run script the same first lines: the data root,
where it came from, and the run's options.

Only the standard library is used, so the flatmap environment can import it
too.
"""

import os
import tomllib

PRODUCTION_CODE = r'D:\sep_histology\code'
PRODUCTION_DATA = r'D:\sep_histology\data'
SNAPSHOT_PREFIX = r'G:\sep_histology_snapshot'


def _canonical(path):
    """Absolute and normalised, with the drive letter in upper case."""
    path = os.path.normpath(os.path.abspath(path))
    drive, rest = os.path.splitdrive(path)
    return drive.upper() + rest


def code_root():
    """The folder that holds get_paths.m, searched upwards from this file."""
    start = os.path.dirname(os.path.abspath(__file__))
    folder = start
    while not os.path.isfile(os.path.join(folder, 'get_paths.m')):
        parent = os.path.dirname(folder)
        if parent == folder:
            raise RuntimeError(f'no get_paths.m in {start} or above it: cannot tell '
                               f'where the code, and so the data, lives')
        folder = parent
    return _canonical(folder)


def _inside(path, folder):
    """Whether `path` is `folder` or lies inside it, whatever the spelling."""
    path, folder = os.path.normcase(path), os.path.normcase(folder)
    return path == folder or path.startswith(folder + os.sep)


def data_root():
    """The data root: SEP_DATA_ROOT if set, else the code folder's sibling data."""
    code = code_root()
    override = os.environ.get('SEP_DATA_ROOT', '')
    if override:
        data = _canonical(override)
    else:
        data = _canonical(os.path.join(os.path.dirname(code), 'data'))
    production_code = os.path.normcase(code) == os.path.normcase(PRODUCTION_CODE)
    if _inside(data, PRODUCTION_DATA) and not production_code:
        raise RuntimeError(f'this copy of the code ({code}) would use the production '
                           f'data ({data}). Set SEP_DATA_ROOT to the data of its own '
                           f'check tree.')
    if os.path.normcase(data).startswith(os.path.normcase(SNAPSHOT_PREFIX)):
        raise RuntimeError(f'the data root {data} is inside the snapshot on G:, which '
                           f'is a backup and never a data root')
    if _inside(data, PRODUCTION_CODE):
        raise RuntimeError(f'the data root {data} is inside the code folder (a copy of the '
                           f'code in a worktree there?). Set SEP_DATA_ROOT to its check tree.')
    return data


DATA = data_root()

# the settings file sits next to the package, beside the run scripts
_SETTINGS_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              'settings.toml')
with open(_SETTINGS_PATH, 'rb') as _fh:
    SETTINGS = tomllib.load(_fh)


def print_settings(options):
    """Print the settings a run works with, as the first lines of its log.

    The data root and where it came from (SEP_DATA_ROOT, or beside the code
    folder), the settings file, then the run's own options, one per line, from
    `options` ({name: value}).
    """
    origin = 'SEP_DATA_ROOT' if os.environ.get('SEP_DATA_ROOT', '') else 'beside the code folder'
    print(f'data root  {DATA}  ({origin})', flush=True)
    print(f'settings   {_SETTINGS_PATH}', flush=True)
    for name, value in options.items():
        print(f'{name:10s} {value}', flush=True)
