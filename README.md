# README

This code project uses Cobaya with cosmological data including SNIa, OHD, and DESI BAO to constrain cosmography parameters via MCMC in support of the paper “Model-independent cosmology for undergraduates: Cosmography, MCMC and the Hubble tension.”

## Project Structure

This project mainly consists of four parts:

```text
.
├── code/      # Program code
├── data/      # Data files
├── chains/    # MCMC chains
└── results/   # Results from running figures.ipynb
```

- `code/`: Stores program code.
- `data/`: Stores data files.
- `chains/`: Stores MCMC chains.
- `results/`: Stores the results obtained by running `figures.ipynb`.

## Description of the `code/` Folder

### `Pantheon+.py`

`code/Pantheon+.py` is used to process the SNIa Pantheon+ sample data files. This script will:

- Remove the SH0ES-calibrated supernova data;
- Apply a redshift cutoff, keeping data with `z > 0.01`.

Both the new and old Pantheon+ data are stored in the `data/` folder.

### `likelihood.py`

`code/likelihood.py` defines the likelihood functions used to constrain the Cosmography parameters with the following observational data:

- SNIa
- OHD
- BAO

### Main Program

The main program is used for MCMC sampling with Cobaya. The main program contains the following three files:

- `main.ipynb`
- `main.py`
- `main.yaml`

You can choose any one of them to run according to your preference.

Before running the main program, Cobaya needs to be installed in advance. For installation instructions, please refer to the official documentation:

https://cobaya.readthedocs.io/en/latest/installation.html

### `figures.ipynb`

`code/figures.ipynb` contains the plotting code, which is used to plot:

- Contour plots;
- The posterior distribution of the Hubble constant.

After running `figures.ipynb`, the generated results are saved in the `results/` folder.

## Suggested Running Workflow

1. Install Cobaya by referring to the official installation documentation.
2. Make sure the required data are prepared in `data/`; if necessary, first run `code/Pantheon+.py` to process the Pantheon+ data.
3. Choose one of `main.ipynb`, `main.py`, and `main.yaml` to run the MCMC sampling.
4. The MCMC chains are stored in `chains/`.
5. Run `code/figures.ipynb` to plot the contour plots and the Hubble constant posterior distribution. The results are saved in `results/`.
