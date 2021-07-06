from compute_model_quants import *

# Length of the sequence
L = 100
# The matrix of effect sizes (betas):
gc_beta = np.repeat([[0., 1., 1., 0.]], L, axis=0).T / L
# The ancestral sequence:
anc_seq = "".join(np.random.choice(list('ACGT'), L))
# The stationary probabilities of the 4 DNA letters:
pi = [.28, .22, .22, .28]

# Compute Zeq:
zeq = compute_z_eq(gc_beta, pi)
print("Zeq:", zeq)

# Compute Z0:
z0 = compute_z0(anc_seq, gc_beta)
print("Z0:", z0)

# Compute sigma_eq:
sigma_eq = compute_sigma_eq(gc_beta, pi)
print("Sigma Eq:", sigma_eq)

# Compute Psi:
psi = compute_psi(anc_seq, gc_beta, pi)

print("Psi (full computation):", psi)
print("Psi (Analytical function of Z0/Zeq and L):", (1./L)*((1. - zeq)*z0 + zeq*(1. - z0)))
