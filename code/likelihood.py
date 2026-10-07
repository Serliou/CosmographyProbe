import numpy as np

# Physical constants
c = 2.99792458e5  # speed of light [km/s]
rd = 147.09        # sound horizon [Mpc], from Planck 2018 [arxiv:1807.06209]

# ----------------------------------------------------------------------
# Cosmography functions
# ----------------------------------------------------------------------
def H(z, H0, q0, j0, s0):
    """ Hubble parameter H(z) [km/s/Mpc] """
    y = z / (1 + z)
    Hoy = H0 * (1 + (1 + q0) * y + (1 + q0 + j0/2. - q0**2/2.) * y**2 +
                (6 + 3*j0 + 6*q0 - 4*q0*j0 - 3*q0**2 + 3*q0**3 - s0) * y**3 / 6.
                + (1 + q0 - 2*j0*q0 + 3*q0**3/2. - s0/2.) * y**4)
    return Hoy

def dL(z, H0, q0, j0, s0):
    """ Luminosity distance [Mpc] """
    y = z / (1 + z)
    dLoy = c / H0 * (y + (3 - q0) * y**2 / 2.
                     + (11 - j0 - 5*q0 + 3*q0**2) * y**3 / 6.
                     + (50 - 7*j0 - 26*q0 + 10*q0*j0 + 21*q0**2 - 15*q0**3 + s0) * y**4 / 24.)
    return dLoy

def dM_over_rd(z, H0, q0, j0, s0):
    """ Transverse comoving distance / sound horizon (for BAO) """
    dM = dL(z, H0, q0, j0, s0) / (1. + z)
    return dM / rd

def dH_over_rd(z, H0, q0, j0, s0):
    """ c / H(z) / rd (for BAO) """
    return c / H(z, H0, q0, j0, s0) / rd

def dV_over_rd(z, H0, q0, j0, s0):
    """ Angle-averaged distance / rd (for BAO) """
    dl = dL(z, H0, q0, j0, s0)
    Hz = H(z, H0, q0, j0, s0)
    dV = (c * z * dl**2 / (1. + z)**2 / Hz) ** (1./3.)
    return dV / rd


# Omit constant term such as ln(C/2\pi), ln(2\pi\sigma^2) in log-likelihood: no impact on final results.

# -------------------------------------------------------------
# SNIa data from the Pantheon+ sample [arxiv:2202.04077]
# -------------------------------------------------------------
# Load Pantheon+ data (with SH0ES calibrators removed)
data = np.genfromtxt('../data/Pantheon+.dat', delimiter=None, names=True,
                     dtype=None, encoding='ascii')
zsn = data['zHD']                 # redshift (CMB frame, recommended by the collaboration)
m_obs = data['m_b_corr']          # observed corrected apparent magnitude
nsn = len(zsn)                    # number of SNe after removing SH0ES calibrators

# Read covariance matrix (first line: dimension, then one element per line)
with open('../data/Pantheon+_STAT+SYS.cov', 'r') as f:
    dim = int(f.readline().strip())
cov_flat = np.loadtxt('../data/Pantheon+_STAT+SYS.cov', skiprows=1)
covsn = cov_flat.reshape(dim, dim)   # full statistical + systematic covariance

# Precompute inverse covariance (use pseudo-inverse for numerical stability)
invcov = np.linalg.pinv(covsn)       # can also use np.linalg.inv if well-conditioned

# SNIa log-likelihood (marginalising over absolute magnitude M_B)
def snlike(H0, q0, j0, s0):
    """
    Compute log-likelihood (ln L) for the Pantheon+ sample.
    Input: H0 [km/s/Mpc], q0, j0, s0 (Cosmography parameters)
    Returns: ln(L) = -0.5 * chi2, where chi2 is analytically marginalised over M_B.
    """
    # Compute theoretical luminosity distance [Mpc]
    dl = dL(zsn, H0, q0, j0, s0)   # actual distance in Mpc
    # Theoretical apparent magnitude (without M_B and constant 25, absorbed by marginalisation)
    th_mag = 5.0 * np.log10(dl)      # i.e., 5 log10(dL/Mpc). Omit constant 5log10(H0/c); absorbed by M' with M, so using dL directly is equivalent.
    delta = m_obs - th_mag           # residual vector

    # Chi-squared after marginalising over M_B
    Acov = np.dot(delta, invcov)          # Δ^T C^{-1}
    magA = np.dot(Acov, delta)            # Δ^T C^{-1} Δ
    magB = np.sum(Acov)                   # Δ^T C^{-1} 1
    magC = np.sum(invcov)                 # 1^T C^{-1} 1
    chi2 = magA - magB**2 / magC

    return -0.5 * chi2


# --------------------------------------------------
# OHD data from Table 1 in [arxiv:2408.02536]
# --------------------------------------------------
# Load OHD data
ccdata = np.loadtxt('../data/OHD.data', unpack=True)
zcc = ccdata[0]
Hcc = ccdata[1]
ercc = ccdata[2]

# OHD log-likelihood
def cclike(H0, q0, j0, s0):   
    H_theory = H(zcc, H0, q0, j0, s0) 
    chi2 = np.sum(((H_theory - Hcc) / ercc) ** 2)
    return -0.5 * chi2


# -------------------------------------------
# DESI-BAO data from [arxiv:2408.02536]
# -------------------------------------------
# Define BAO data: (redshift, distance type, observed value, error)
bao_data = [
    (0.30, 'dV', 7.93, 0.15),
    (0.51, 'dM', 13.62, 0.25),
    (0.51, 'dH', 20.98, 0.61),
    (0.71, 'dM', 16.85, 0.32),
    (0.71, 'dH', 20.08, 0.60),
    (0.93, 'dM', 21.71, 0.28),
    (0.93, 'dH', 17.88, 0.35),
    (1.32, 'dM', 27.79, 0.69),
    (1.32, 'dH', 13.82, 0.42),
    (1.49, 'dV', 26.07, 0.67),
    (2.33, 'dM', 39.71, 0.94),
    (2.33, 'dH', 8.52, 0.17),
]
# Convert to arrays for indexing
z_bao = np.array([d[0] for d in bao_data])
types = [d[1] for d in bao_data]
obs = np.array([d[2] for d in bao_data])
err = np.array([d[3] for d in bao_data])

# BAO log-likelihood
def baolike(H0, q0, j0, s0):   
    # Pre-allocate theoretical values
    theory = np.zeros_like(obs)
    for i, (z, typ) in enumerate(zip(z_bao, types)):
        if typ == 'dV':
            theory[i] = dV_over_rd(z, H0, q0, j0, s0)
        elif typ == 'dM':
            theory[i] = dM_over_rd(z, H0, q0, j0, s0)
        elif typ == 'dH':
            theory[i] = dH_over_rd(z, H0, q0, j0, s0)
    chi2 = np.sum(((theory - obs) / err) ** 2)
    return -0.5 * chi2
