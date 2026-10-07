"""
Extract Pantheon+ sample without SH0ES calibrators and apply a redshift cut.

This script reads the official Pantheon+SH0ES data release files:
- Data: Pantheon+SH0ES.dat (with header)
- Covariance: Pantheon+SH0ES_STAT+SYS.cov (first line: dimension, then elements)

It selects only SNe with:
    IS_CALIBRATOR == 0  (non-calibrator Hubble-flow SNe)
    zHD > 0.01          (remove low-redshift SNe dominated by peculiar velocities)

It saves the filtered data with its original header as 'Pantheon+.dat',
and writes the corresponding sub-covariance matrix with its dimension header
as 'Pantheon+_STAT+SYS.cov'.

Usage:
    python extract_pantheon.py
"""

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------
# File paths (modify as needed)
DATA_FILE = '../data/Pantheon+SH0ES.dat'
COV_FILE = '../data/Pantheon+SH0ES_STAT+SYS.cov'
OUT_DATA = '../data/Pantheon+.dat'
OUT_COV = '../data/Pantheon+_STAT+SYS.cov'

# ----------------------------------------------------------------------
def main():
    # 1. Load the data file (with header)
    print(f"Reading data from {DATA_FILE} ...")
    data = pd.read_csv(DATA_FILE, delim_whitespace=True, header=0, comment='#')
    print(f"Total number of SNe: {len(data)}")
    print(f"Columns: {data.columns.tolist()}")

    # 2. Load the covariance matrix
    print(f"Reading covariance matrix from {COV_FILE} ...")
    with open(COV_FILE, 'r') as f:
        lines = f.readlines()
    # First line is the dimension
    dim = int(lines[0].strip())
    print(f"Original covariance dimension: {dim}")
    # Remaining lines are matrix elements (one per line)
    elements = [float(line.strip()) for line in lines[1:] if line.strip()]
    # Check length
    if len(elements) != dim * dim:
        raise ValueError(f"Expected {dim*dim} elements, got {len(elements)}")
    cov_full = np.array(elements).reshape(dim, dim)
    print(f"Covariance matrix shape: {cov_full.shape}")

    # 3. Create mask for:
    #    - non-calibrator SNe (IS_CALIBRATOR == 0)
    #    - redshift cut zHD > 0.01
    mask = (data['IS_CALIBRATOR'] == 0) & (data['zHD'] > 0.01)
    selected_data = data[mask].copy()
    selected_indices = np.where(mask)[0]   # row indices in the original data

    print(f"Number of selected SNe (non-calibrator, zHD > 0.01): {len(selected_data)}")

    # 4. Extract the sub-covariance matrix
    cov_sub = cov_full[np.ix_(selected_indices, selected_indices)]
    print(f"Sub-covariance matrix shape: {cov_sub.shape}")

    # 5. Save the filtered data with header (same format as input)
    print(f"Writing selected data to {OUT_DATA} ...")
    selected_data.to_csv(OUT_DATA, sep=' ', index=False, header=True, float_format='%.8f')
    print(f"Data saved with header.")

    # 6. Save the sub-covariance matrix with dimension header
    print(f"Writing sub-covariance to {OUT_COV} ...")
    new_dim = cov_sub.shape[0]
    with open(OUT_COV, 'w') as f:
        f.write(f"{new_dim}\n")
        # Write each element on a new line, in row-major order
        for val in cov_sub.ravel():
            f.write(f"{val:.8f}\n")
    print(f"Covariance saved with new dimension header.")

    print("Done.")

if __name__ == '__main__':
    main()
