# TTL_pipeline

A metabolomics data-processing pipeline, currently built as a set of
standalone Python modules — one per pipeline stage — with a shared
environment during this early development phase. The plan is to eventually
wire these modules together as a Nextflow/nf-core pipeline with a GUI front
end; for now each module is a plain, independently testable Python script.

## Modules

| Module | Status | Description |
|---|---|---|
| [`modules/Upload_and_validator_module`](modules/Upload_and_validator_module/README.md) | Done | Validates mzML files, sample metadata, and a data matrix (required), plus optional feature table / MGF / annotation list, before downstream processing. |

New modules are added as `modules/<module_name>/`, each with its own
script, `tests/`, and a module-level `README.md` explaining what it does
and how to run it.

## Setup

All modules currently share a single environment. If you're using conda:

```bash
conda create -n TTL_pipeline python=3.11 pandas pytest -y
conda activate TTL_pipeline
pip install -r requirements.txt
```

Or with a plain virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate        # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

`requirements.txt` at the project root covers dependencies for every
module; as modules are added, their dependencies get appended here rather
than each module carrying its own `requirements.txt`. If a module ever
needs a conflicting dependency version, that's the point at which it gets
split into its own environment — not before.

## Running a module

Each module is run from inside its own folder — see that module's README
for exact usage. For example:

```bash
cd modules/upload_validation
./run_check.sh
```

## Running tests

Each module has its own `tests/` folder, runnable independently:

```bash
cd modules/upload_validation
pytest -v
```
