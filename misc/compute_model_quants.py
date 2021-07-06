"""
Author: Shadi Zabad
July 2021

This script implements functions for estimating the central quantities in the
Neutral Model outlined in Zabad & Moses (2021).
"""

import numpy as np
import itertools


def compute_psi(ancestral_seq, beta, pi, alphabet=('A', 'C', 'G', 'T')):
    """
    This function takes three arguments:
    - ancestral_seq: DNA sequence (length L) representing the ancestral sequence for a set of species
    - beta: A 4xL matrix of effect sizes for each position/allele in the sequence
    - pi: The stationary distribution of the alleles

    NOTE: If your alleles are not DNA-based, make sure to provide an alphabet as well.

    And returns an estimate for the parameter Psi.
    """

    # The index of each fixed allele in the ancestral sequence:
    idx = [alphabet.index(l) for i, l in enumerate(ancestral_seq)]

    psi = 0.

    for j in range(len(ancestral_seq)):
        for k in range(4):
            if k != idx[j]:
                psi += pi[k] * (beta[idx[j], j] - beta[k, j]) ** 2

    return psi


def compute_sigma_eq(beta, pi):
    """
    This function takes two arguments:
    - beta: A 4xL matrix of effect sizes for each position/allele in the sequence
    - pi: The stationary distribution of the DNA letters (ACGT, in alphbetical order)

    And returns an estimate for the parameter Psi.
    """

    sig_eq = 0.

    for j in range(beta.shape[1]):
        for k, m in itertools.combinations(range(len(pi)), 2):
            sig_eq += pi[k] * pi[m] * (beta[k, j] - beta[m, j]) ** 2

    return sig_eq


def compute_z_eq(beta, pi):
    """
    This function takes two arguments:
    - beta: A 4xL matrix of effect sizes for each position/allele in the sequence
    - pi: The stationary distribution of the alleles

    And returns an estimate of the stationary value for the mean phenotype.
    """

    pi_mat = np.repeat([pi], beta.shape[1], axis=0).T

    return (beta * pi_mat).sum()


def compute_z0(ancestral_seq, beta, alphabet=('A', 'C', 'G', 'T')):
    """
    This function takes two arguments:
    - ancestral_seq: DNA sequence (length L) representing the ancestral sequence for a set of species
    - beta: A 4xL matrix of effect sizes for each position/allele in the sequence

    NOTE: If your alleles are not DNA-based, make sure to provide an alphabet as well.

    And returns an estimate of the ancestral mean phenotype.
    """

    # The index of each fixed allele in the ancestral sequence:
    idx = [alphabet.index(l) for i, l in enumerate(ancestral_seq)]

    # The ancestral sequence in one-hot encoding:
    one_hot = np.zeros(shape=(4, len(ancestral_seq)))
    one_hot[idx, np.arange(len(ancestral_seq))] = 1

    return (one_hot * beta).sum()

