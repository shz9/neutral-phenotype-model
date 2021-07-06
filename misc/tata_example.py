from compute_model_quants import *

# Length of the sequence
L = 10

# The transcription factor binding probability matrix for the TATA box:
tata_tf_bpm = np.array([
    [0.157760814, 0.371501272, 0.389312977, 0.081424936],
    [0.043256997, 0.119592875, 0.048346056, 0.788804071],
    [0.89821883, 0.002544529, 0.007633588, 0.091603053],
    [0.010178117, 0.027989822, 0.007633588, 0.954198473],
    [0.903307888, 0.002544529, 0.015267176, 0.078880407],
    [0.684478372, 0.002544529, 0.002544529, 0.31043257],
    [0.942558747, 0.010443864, 0.028720627, 0.018276762],
    [0.567430025, 0.007633588, 0.114503817, 0.31043257],
    [0.396946565, 0.114503817, 0.402035623, 0.086513995],
    [0.145038168, 0.34605598, 0.384223919, 0.124681934]
])

# The stationary probabilities of the 4 DNA letters:
pi = np.array([0.32, 0.18, 0.18, 0.32])

tata_beta = np.log(tata_tf_bpm / pi).T
# The ancestral sequence:
anc_seq = "".join(np.random.choice(list('ACGT'), L))

# Compute Zeq:
zeq = compute_z_eq(tata_beta, pi)
print("Zeq:", zeq)

# Compute Z0:
z0 = compute_z0(anc_seq, tata_beta)
print("Z0:", z0)

# Compute sigma_eq:
sigma_eq = compute_sigma_eq(tata_beta, pi)
print("Sigma Eq:", sigma_eq)

# Compute Psi:
psi = compute_psi(anc_seq, tata_beta, pi)
print("Psi:", psi)
