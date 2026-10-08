# polySIM

Build hydrogen-terminated polymer chains from two-ended pSMILES, generate an
MMFF conformer, optionally refine and analyze the isolated chain with xTB, and
pack rigid chain copies into an approximate periodic simulation box.

> **Project status:** research prototype. The packed multi-chain box is a
> starting structure, not an equilibrated molecular-dynamics system.

This root-level package converts a two-ended repeat-unit pSMILES into:

1. a hydrogen-terminated finite chain;
2. an RDKit MMFF94/MMFF94s conformer;
3. optionally, a GFN-FF geometry-optimized and xTB-analyzed **single chain**; and
4. an unoptimized approximate periodic box containing rigid chain copies.

Linear-chain assembly uses `stk.polymer.Linear`, following the supplied STK
building-block approach. RDKit remains responsible for pSMILES parsing,
STK/RDKit interchange, and MMFF conformer generation. xTB is never applied to
the assembled multi-chain system.

For `N >= 2`, STK receives an H/Br initial unit (`A`), `N-2` Br/Br middle
units (`B`), and a Br/H terminal unit (`C`). Thus the topology sequence is
`A + B*(N-2) + C`, containing exactly `N` chemical repeat units. For `N=1`,
both wildcard sites are replaced by hydrogen directly.

## Requirements

PolySIM is an installable Python package. It requires Python 3.11 or newer.
Python 3.11 and 3.12 are the intended runtime versions; the tested STK release
requires Python 3.11+.

The simplest complete setup uses Conda and includes the external xTB program:

```bash
conda env create -f environment.yml
conda activate polisim
python -m pytest
```

Alternatively, create a virtual environment and install the Python dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate             # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Geometry optimization also requires the external `xtb` executable. It is not a
Python dependency and is therefore not installed by `pip`; one supported route
is:

```bash
conda install -c conda-forge xtb
xtb --version
```

The full workflow is tested with STK 2026.1.4.0 and xTB 6.5.1. The chain and
MMFF stages can run without xTB by passing `--xtb-method none`.

For development and testing, install the `dev` dependency group:

```bash
python -m pip install -e ".[dev]"
python -m pytest
```

## Quick start

Run a small workflow that does not require xTB:

```bash
python -m polysim \
  --psmiles-file examples/pvdf.psmiles.txt \
  --repeat-units 5 \
  --chains 3 \
  --xtb-method none \
  --output-directory output/quickstart
```

Inspect every available option with `python -m polysim --help`. After
installation, the equivalent shorter command is `polysim --help`.

## pSMILES input

Place exactly one pSMILES on a non-empty, non-comment line. Explicit head and
tail positions are determined by their order in the string:

```text
# pvdf.psmiles.txt
[*]CC(F)(F)[*]
```

An example is provided at `polySIM/examples/pvdf.psmiles.txt`.

The builder does not canonicalize the pSMILES before interpreting its connection
points. RDKit input atom order is retained: the first `[*]` is the head and the
second is the tail. After exactly `N` repeat units are connected, the two outer
valences are hydrogen-capped. For example, `N=10` PVDF gives `C20H22F20`.

By default, the two wildcard bonds must have the same bond order. Some PI1M
records use different head/tail bond orders. For those, explicitly choose
`--mismatched-bond-policy single` to reproduce the supplied STK builder's
effective single-bond linking behavior. The default remains `strict` so
bond-order changes are never silent.

## Binary random copolymers

Supply a second two-ended monomer and its percentage to build one binary random
copolymer. `--repeat-units` is always the total number of units in the chain and
may differ between runs. For example, this creates a 17-unit chain containing
approximately 30% monomer 2:

```bash
python -m polysim \
  --monomer-1-file examples/monomer_1.psmiles.txt \
  --monomer-2-file examples/monomer_2.psmiles.txt \
  --monomer-2-percentage 30 \
  --repeat-units 17 \
  --chains 10 \
  --random-seed 42 \
  --xtb-method none \
  --output-directory output/random_copolymer
```

PolySIM calculates the nearest integer number of monomer-2 units using
round-half-up, creates the exact required counts, and then shuffles them. Thus a
17-unit chain requested at 30% contains 5 monomer-2 units and realizes 29.41%.
The same seed and inputs reproduce the same sequence. Values of 0% and 100%
produce pure monomer-1 and pure monomer-2 chains, respectively.

The requested percentage, realized percentage, counts, seed, and complete A/B
sequence are saved in `manifest.json`; the sequence is also written to
`chain_sequence.txt`. The current packing stage uses rigid copies of this one
realized chain, so all chains in `system.xyz` have the same sequence and length.
It does not yet create a different sequence or a different repeat-unit count for
each packed chain.

## One-command workflow

Run the workflow from inside the `polySIM` folder:

```bash
cd polySIM
python -m polysim \
  --psmiles-file examples/pvdf.psmiles.txt \
  --repeat-units 50 \
  --chains 30 \
  --xtb-method gfn2 \
  --density 0.3 \
  --minimum-distance 2.0 \
  --output-directory output/pvdf_30x50
```

The default is a staged single-chain calculation: GFN-FF geometry optimization
with xTB's inertial optimizer, followed by a GFN2-xTB single-point energy on the
relaxed geometry. Charge 0, UHF 0, no solvent, 14 OpenMP threads, and no Hessian
follow the supplied reference scripts. The inertial optimizer is used because
the default approximate-normal-coordinate optimizers can fail while constructing
a model Hessian for long, flexible chains.

Override the stages with `--xtb-method`, `--xtb-geometry-method`,
`--xtb-optimizer-engine`, `--xtb-charge`, `--xtb-uhf`, `--xtb-threads`,
`--xtb-opt-level`, and `--xtb-hessian`. The Hessian option applies to the
geometry-optimization stage.

Supported `--xtb-method` values:

- `gfn2` (default), `gfn1`, or `gfn0` for the post-optimization electronic
  single point;
- `gfnff` to run only the GFN-FF geometry optimization;
- `none` to pack the converged MMFF chain without xTB.

Set `--xtb-geometry-method none` to omit the GFN-FF relaxation and run the
selected electronic method as a single point on the MMFF geometry. xTB is used
only for one isolated chain. A failure stops the workflow and retains the
`.xtb.log`; it does not silently substitute a result.

## xTB analysis

After xTB, the workflow writes `xtb_analysis/xtb_analysis.json` and `.csv`.
The analysis extracts the final total energy, reports relative energies, ranks
results, and stages selected geometries. The xTB analysis module can also
analyze a directory containing several same-atom-order xTB conformers. It uses ordered
Kabsch RMSD filtering, avoiding unreliable bond inference from XYZ files.

Standalone usage:

```bash
python -m polysim.xtb.analysis \
  --input-directory path/to/xtb/results \
  --output-directory path/to/analysis \
  --rmsd-threshold 2.0 \
  --max-selected 3
```

## Output

- `input.psmiles.txt`: copied input for provenance;
- `input.monomer_2.psmiles.txt`: copied second input for copolymer runs;
- `chain_sequence.txt`: realized A/B sequence for copolymer runs;
- `chain_topology.sdf`: finite-chain topology;
- `chain_mmff.xyz`: MMFF-optimized chain;
- `chain_mmff_explicit_h.sdf`: explicit-H topology matching XYZ atom order;
- `chain_gfnff_opt.xyz` and `.xtb.log` for the geometry optimization;
- `chain_<method>_sp.xyz` and `.xtb.log` for an electronic single point;
- `xtb_analysis/`: xTB energy and conformer-selection reports;
- `system.xyz`: all chains in one approximate box;
- `system.box.json`: box dimensions and atom ranges for every chain;
- `manifest.json`: inputs, methods, convergence, and output inventory.

XYZ does not contain bonds or periodic box vectors. For handover, keep
`system.xyz`, `system.box.json`, `chain_mmff_explicit_h.sdf`, and
`manifest.json` together.

## Code organization

```text
polysim/chemistry/  -> RDKit/STK chain construction and MMFF conformers
polysim/xtb/        -> isolated-chain xTB execution and result analysis
polysim/md/         -> approximate box construction for downstream MD
polysim/io/         -> shared molecular file readers and writers
polysim/workflow.py -> complete orchestration and command-line interface
```

Each section can be inspected and tested independently. The primary supported
interface is `python -m polysim`; specialist module commands include
`python -m polysim.xtb.analysis --help` and
`python -m polysim.xtb.calculator --help`.

## Repository layout

```text
polySIM/
├── polysim/
│   ├── chemistry/
│   │   ├── chains.py        # pSMILES and STK chain construction
│   │   └── conformers.py    # RDKit MMFF conformer generation
│   ├── xtb/
│   │   ├── calculator.py    # xTB execution wrapper
│   │   └── analysis.py      # energy and RMSD analysis
│   ├── md/
│   │   └── system.py        # rigid-chain packing / MD preparation
│   ├── io/
│   │   └── xyz.py           # shared XYZ reader and writer
│   ├── workflow.py          # end-to-end workflow
│   └── __main__.py          # python -m polysim entry point
├── examples/                # example pSMILES inputs
├── tests/                   # fast workflow tests
├── pyproject.toml           # package metadata and all Python dependencies
└── environment.yml          # complete Conda environment, including xTB
```

## Testing

The test suite exercises repeat-unit assembly, hydrogen termination, random
copolymer composition, MMFF conformer generation, packing constraints, and xTB
output parsing. It does not invoke the external xTB executable.

```bash
python -m pytest
```

## Reproducibility and limitations

- Random sequence generation, conformer generation, and packing use
  `--random-seed` (default: `42`). The manifest records workflow settings.
- All packed chains are rigid copies of one generated chain and, for random
  copolymers, share the same sequence.
- The packer creates an approximate cubic box at the requested density. It does
  not minimize or equilibrate the assembled system.
- XYZ files do not retain bonds or periodic box vectors; keep the accompanying
  SDF, box metadata, and manifest files.
- Calculation outputs are ignored by Git by default. Commit selected results
  explicitly with `git add -f <path>` only when they are intentional fixtures.

## Contributing

Create a focused branch, add or update tests for behavioral changes, run
`python -m pytest`, and keep generated calculation data out of commits. Please
open an issue before making changes that alter chemical assumptions or file
formats.

## License

No license has been selected yet. Until a `LICENSE` file is added, the code is
not granted for redistribution or reuse outside the rights provided by default
copyright law.
