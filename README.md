# SEP-GluA1 whole-brain maps: analysis code

![Slice-order montage of a P16 brain](assets/slice_order_montage.png)

*A P16 brain (MG911) in the slice-order editor of `run_order_slices`, six of
its eight columns: each section labelled with its place in the curated order,
its index as extracted in brackets, F if flipped; red = DAPI, green =
nanobody, blue = autofluorescence.*

Analysis code for whole-brain maps of surface GluA1, an AMPA receptor
subunit, in SEP-GluA1 knock-in mice. Coronal sections are stained without
detergent, or with very little, with a GFP-booster nanobody that binds the
SEP tag, so in principle only receptors at the surface are labelled; no
control has tested this on these brains yet. Each section is imaged in four
channels (DAPI, the nanobody, autofluorescence, and the green of SEP itself)
and registered to an atlas of the mouse's age: the Allen CCF for adults, the
DeMBA atlas of their age (Carey 2025) for young mice.

The project serves Aim 2.1 of the SNSF Weave grant *Dendritic plasticity
rules shaping the emergence of multisensory integration in the developing
cortex*: surface GluA1 as a proxy for synaptic plasticity potential across
development. The dendritic imaging of the same grant is in the two-photon
imaging repository (`D:\dendrites\code`). See
[`docs/SCIENTIFIC_CONTEXT.md`](docs/SCIENTIFIC_CONTEXT.md) for the questions,
the results and the papers, and [`docs/ROADMAP.md`](docs/ROADMAP.md) for
what is open.

## The four lines of work

1. **Where an experience changes the map.** Naive mice against mice after
   rhythmic whisker stimulation or after a detection task. Paused, not
   closed. MATLAB, `group_comparison/`.
2. **How the signal is distributed in the adult brain**, and how
   reproducibly across mice. Python, `mapping/`; the MATLAB `P8` and `P10`
   until the Python route answers their questions.
3. **What the map measures**, against the Allen in situ hybridisation maps
   and the green SEP channel. Python, `mapping/`; `P9` until replaced.
4. **How young brains differ from adult ones**, for the grant. Python,
   `mapping/`.

Preprocessing and registration (MATLAB, with LightSuite and an automatic
annotation of control points) serve all four.

## Setup

- **MATLAB** R2024b on Windows (R2022b at least), with the Image Processing,
  Computer Vision, Optimization, Statistics and Machine Learning, and
  Parallel Computing toolboxes. Once per session, in a fresh MATLAB:
  `restoredefaultpath; cd('D:\sep_histology\code'); sep_setup_paths`.
- **Paths.** The code and data folders sit side by side, `<root>\code` and
  `<root>\data` (here `D:\sep_histology\`). `get_paths.m` and
  `mapping/sepmap/config.py` work the data root out from where the code
  sits. The environment variable `SEP_DATA_ROOT` moves the whole data tree,
  inputs and outputs, for a check on a copy; a copy of the code is refused
  the production data.
- **elastix 5.1.0** on the PATH (LightSuite calls it): `elastix --version`
  must answer in a new terminal before MATLAB starts.
- **Python**: Anaconda's Python 3.12, one environment per job, each made
  from a pinned file whose top lines give the commands:

  | environment | runs | made from |
  |---|---|---|
  | `tools\venv_atlas` | every `mapping\run_*.py` but `run_closeup`; `atlas\build_demba_atlas.py`; the Python check tools | `tools\requirements_atlas.txt` |
  | `tools\venv_flat` | `mapping\run_closeup.py`, the cortical flatmaps | `tools\requirements_flat.txt` |
  | `registration\auto_annotation\.venv` | the automatic annotation's engine | `registration\auto_annotation\setup.ps1`; versions pinned in `tools\requirements_auto_annotation.txt` |
  | `tools\venv_dev` | pytest and ruff | `tools\requirements_dev.txt` |

  The engine's `setup.ps1`, run once per machine, installs torch and ends
  with a self-test; without it the control-point GUI has no automatic keys.
- **Atlases**, under the data root: `atlas\` (Allen CCFv3, 10 µm),
  `atlas_demba_p<age>\` (built by `atlas\build_demba_atlas.py <age>`),
  `atlas_flatmap\` and `atlas_ish\`; see [`atlas/README.md`](atlas/README.md).

## Pipelines

| folder | what | start with |
|---|---|---|
| `preprocessing/` | raw `.czi` files to centred, ordered, corrected sections (MATLAB) | [`preprocessing/README.md`](preprocessing/README.md) |
| `registration/` | sections to the atlas of the mouse's age, control points by hand or proposed and reviewed; the SEP channel (MATLAB, Python engine) | [`registration/README.md`](registration/README.md) |
| `group_comparison/` | line 1: naive against RWS or behaviour (MATLAB) | [`group_comparison/README.md`](group_comparison/README.md) |
| `mapping/` | lines 2 to 4: per-brain volumes, the adult map, the ISH comparison, young against adult (Python) | [`mapping/README.md`](mapping/README.md) |
| root: `P8_*.m`, `P9_*.m`, `P10_*.m` | the earlier adult, ISH and autofluorescence analyses (MATLAB) | the header of each script |

A new brain end to end, with the manual steps and the traps:
[`docs/ADDING_DATA.md`](docs/ADDING_DATA.md). Each driver (`run_*.m`,
`mapping\run_*.py`) lists its pipeline's run order in its header. MATLAB
settings sit under `%% Settings` in each driver, the Python route's in
`mapping/settings.toml`; long MATLAB stages run detached, with a log, through
`tools\run_matlab_detached.ps1`.

Shared code:

- `get_paths.m` and `sep_setup_paths.m`: the data root, the MATLAB path
- `common/`: the cohort table (`cohort.csv`), volume reading, left-right
  statistics, the palette; `atlas/`: `get_atlas`, the DeMBA builder, checks
- `tests/`, `tools/` (output and code-identity checks, the detached runner,
  the Python requirements) and `third_party/` (LightSuite with our patches,
  and three smaller packages), each with its README
- `archive/`: retired code, kept until checked; `docs/`: the documents;
  `assets/`: images for this page

## Where outputs go

Under the data root; none is tracked in git.

- `<group>\<mouse>\` (groups `naive`, `rws`, `behavior`, `young`): the
  copied `.czi` files and `lightsuite\`, every stage up to the registered
  volumes.
- `<group>\` and `comparisons\`: the plasticity comparison (the approved
  run: `comparisons\naive_vs_rws\`, `naive_vs_behavior\`), and P8 to P10.
- `comparisons_v2\`: per-brain and cohort volumes, young against adult, the
  diagnostic sheets, and a README on reading their numbers.
- `adult_v2\`: the beyond-abundance analysis, the channels, the ISH tests.

## Status

The code was reorganised in September and October 2026
([`docs/REFACTOR_PLAN.md`](docs/REFACTOR_PLAN.md)); each step was checked
against the code before it on a reference set of brains. Next: six young
brains of P28 to P36, and one declared set of structures for the adult and
ISH analyses (A1 to A5), after which `P8` to `P10` retire. Before reading a
result, see the limits in [`docs/ROADMAP.md`](docs/ROADMAP.md): above all,
no reading is an absolute level, since young and adult brains were imaged in
different sessions, and the surface claim rests on the staining protocol
alone.

## Conventions

- Code style: [`docs/STYLE.md`](docs/STYLE.md). Raw data is only read,
  outputs go under the data root, and a step overwrites only its own outputs.
- A mouse's age is the one in the name of its raw folder: `MG904_SepGluA_P22`
  is P22. Adult folders end in `_Gria1`; the adults were well above P60,
  their exact ages not recorded.
- `P<n>` in an age, a cohort tag or an atlas key is postnatal day n
  (`young_P20`, `demba_p20`); as script names, P8, P9 and P10 are the last
  of the old pipeline steps.
- `main` is the working branch. Tag `grant-2026-09` is the code behind the
  grant figures, `refactor-start` the code before the reorganisation.
