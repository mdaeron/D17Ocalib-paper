# Data and code for "Triple oxygen isotope equilibrium between calcite and water: new constraints over 8–200 °C and their metrological implications" by Daëron et al.

This is the full dataset and codebase for the following manuscript (currently under review): _Triple oxygen isotope equilibrium between calcite and water: new constraints over 8–200 °C and their metrological implications_, by M. Daëron, T. Clauzel, J. Chaillot, M. Pesnin, J. J. Jautzy, B. R. Fosu, and S. Kassi.

## Contents:

* `code`: Python scripts and libraries used to process the data
* `inputs`: "raw" input data
  + `ftir`: FTIR spectra for our synthetic samples, along with calcite and aragonite reference materials
  + `rawdata.csv`: "reduced" VCOF-CRDS measurements (before isotopic normalization) compiled into a single file
  + `runs`:  "raw" VCOF-CRDS measurements used to compile `rawdata.csv`
  + `temperatures.csv`: crystallization temperatures of our samples
* `justfile`: a `just` configuration file (see https://just.systems/man/en) containing all data processing recipes
* `manuscript`: manuscript files, including Typst source and data exports to be used in the text.
* `output`: all numerical and pdf outputs
* `pyproject.toml`: a standardized Python project configuration file
* `readme.md`: this file
* `uv.lock`: a machine-readable specification of all dependencies required to run this code reproducibly

## Usage

* Install [`just`](https://just.systems/man/en/installation.html)
* Install [`uv`](https://docs.astral.sh/uv/getting-started/installation)
* Running `just` from the root of this repository will fully reprocess the whole dataset.
* Partial reprocessing is possible using the recipes defined in `justfile`, e.g., by running `just harvest` to recompile the raw VCOF-CRDS data into `inputs/rawdata.csv`.
