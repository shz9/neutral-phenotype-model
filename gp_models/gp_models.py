from scipy.linalg import cholesky, cho_solve
import scipy.optimize as optim
import numpy as np
import pandas as pd


def get_time_matrices(tree, tips=None):
    """
    This function takes a Biopython tree and returns the
    shared distance matrix (time to most recent common ancestor - MRCA)
    and
    the pairwise distance matrix (pairwise distance since divergence from MRCA)
    """

    tips = tree.get_terminals() if tips is None else tips

    n_tips = len(tips)  # Number of terminal species

    sdist_matrix = np.zeros((n_tips, n_tips))  # Shared distance matrix
    pdist_matrix = np.zeros((n_tips, n_tips))  # Pairwise distance matrix

    for i in range(n_tips):
        for j in range(i, n_tips):
            if i == j:

                pdist_matrix[i, j] = 0.
                sdist_matrix[i, j] = tree.distance(tree.root, tips[i])

            else:

                mrca = tree.common_ancestor(tips[i], tips[j])

                sdist_matrix[i, j] = sdist_matrix[j, i] = tree.distance(tree.root, mrca)
                pdist_matrix[i, j] = pdist_matrix[j, i] = tree.distance(tips[i], tips[j])

    return sdist_matrix, pdist_matrix


class GaussianProcessModel(object):

    def __init__(self, data, tree):

        self.tree = tree
        self.data = data

        tree_terminals = [ts.name for ts in self.tree.get_terminals()]

        self.tips = list(set(data.index).intersection(set([ts.name for ts in self.tree.get_terminals()])))

        self.training_mask = np.arange(len(self.tips))  # use all data points for training by default
        self.test_mask = None

        if len(self.tips) < 2:
            raise Exception("Less than two species were matched \
            between the phylogenetic tree and the data matrix. Ensure that the names are consistent \
            between the two data sources.")

        for t in tree_terminals:
            if t not in self.tips:
                self.tree.prune(t)

        self.data = self.data[self.tips]
        self.sd, self.pd = get_time_matrices(tree, self.tips)

        self.k = None  # Number of parameters of the model
        self.m = len(self.tips)  # Number of species
        self.n = self.m  # Total number of measured characters
        self.Ky = None

    def mean(self, test=False):
        raise NotImplementedError

    def cov(self, test=False):
        raise NotImplementedError

    def fit_loco(self):

        mse = 0.

        for t in self.tips:
            self.training_mask = [i for i, sp in enumerate(self.tips) if sp != t]
            self.fit()
            self.predict(t)
            mse += (self.predict(t).values - self.data[t])**2

        return {
            'MSE': mse,
            'RMSE': np.sqrt(mse)
        }

    def fit(self):
        raise NotImplementedError

    def predict(self, species=None):

        if species is not None:
            self.test_mask = [i for i, sp in enumerate(self.tips) if sp in species]
        else:
            self.test_mask = self.training_mask

        return pd.Series(self.mean(test=True) + np.dot(self.cov(test=True), self.Ky),
                         index=np.array(self.tips)[self.test_mask])

    def nll(self):
        # The negative loglikelihood

        # Compute the difference between the mean phenotype for each
        # species (x_i) and the expected value under the model (mu_i):
        # (x_i - mu_i)

        mean_diff = self.data[np.array(self.tips)[self.training_mask]] - self.mean()

        # Compute the cholesky factor:
        try:
            L = cholesky(self.cov(), lower=True)
        except Exception as e:
            raise e

        # Compute the log-determinant using the cholesky factor:
        log_det = 2 * np.log(np.diag(L)).sum()

        # Compute the product of the squared mean difference and the precision matrix:
        # (x - mu)'Phi(x - mu)
        # where x is a vector of mean phenotype for each species
        # and mu is the expected value for each species
        # and Phi is the precision matrix (inverse of covariance)

        self.Ky = cho_solve((L, True), mean_diff)
        prod = np.dot(mean_diff, self.Ky)

        return .5 * (
                self.m * np.log(2. * np.pi) + log_det + prod
        )

    def aic(self, nll=None):
        # Akaike information criterion
        if nll is None:
            return 2. * (self.nll() + self.k)
        else:
            return 2. * (nll + self.k)

    def corrected_aic(self, nll=None):
        # Akaike information criterion (corrected for sample size)
        return self.aic(nll) + 2. * (self.k * (self.k + 1.)) / (self.n - self.k - 1.)

    def bic(self, nll=None):
        # Bayesian information criterion
        if nll is None:
            return 2. * self.nll() + self.k * np.log(self.n)
        else:
            return 2. * nll + self.k * np.log(self.n)


class NeutralModel(GaussianProcessModel):

    def __init__(self, data, tree, u=None):

        super().__init__(data, tree)

        self.fixed_u = u is not None

        # The 5 parameters of the model:
        self.z0 = None  # The ancestral phenotype
        self.z_eq = None  # The equilibrium value of the phenotype
        self.psi = None  # Psi, captures the relationship between ancestral and equilibrium sequence
        self.sigma_eq = None  # The equilibrium variance
        self.u = u  # The scaling factor u, or the rate in our context

        self.k = 4 if self.fixed_u else 5  # Number of parameters under the model

    def mean(self, test=False):
        if test:
            sd = np.diag(self.sd)[self.test_mask]
        else:
            sd = np.diag(self.sd)[self.training_mask]

        return self.z0 * np.exp(-self.u * sd) + (1. - np.exp(-self.u * sd)) * self.z_eq

    def cov(self, test=False):
        if test:
            row_mask = self.test_mask
        else:
            row_mask = self.training_mask

        sd = self.sd[np.ix_(row_mask, self.training_mask)]
        pd = self.pd[np.ix_(row_mask, self.training_mask)]
        return np.exp(-self.u * pd) * (1. - np.exp(-self.u * sd)) * (
                np.exp(-self.u * sd) * (self.psi - self.sigma_eq) + self.sigma_eq
        )

    def fit(self):
        # Fit the hyperparameters to the data

        init_params = np.random.uniform(size=self.k)

        bounds = [
            (None, None),
            (None, None),
            (1e-12, None),
            (1e-12, None)
        ]

        if not self.fixed_u:
            bounds += [(1., 2.)]

        def objective(params):
            if self.k == 4:
                self.z0, self.z_eq, self.psi, self.sigma_eq = params
            else:
                self.z0, self.z_eq, self.psi, self.sigma_eq, self.u = params

            return self.nll()

        optim_res = optim.minimize(objective,
                                   init_params,
                                   bounds=bounds)

        return {
            'Optimization': {
                'Success': optim_res.success,
                'Message': optim_res.message
            },
            'Loglikelihood': -optim_res.fun,
            'DOF': self.k,
            'AIC': self.aic(optim_res.fun),
            'AIC.c': self.corrected_aic(optim_res.fun),
            'BIC': self.bic(optim_res.fun),
            'Parameters': dict(zip(['Z0', 'Zeq', 'Psi', 'sigma_eq', 'u'],
                                   optim_res.x))
        }


class OU(GaussianProcessModel):

    def __init__(self, data, tree, equilibrium_z0=False):

        super().__init__(data, tree)

        self.eq_z0 = equilibrium_z0  # set the ancestral phenotype to the equilibrium phenotype

        # The 4 parameters of the model:
        self.z0 = None  # The ancestral phenotype
        self.z_eq = None  # The equilibrium value of the phenotype
        self.sigma = None  # The variance
        self.alpha = None  # The strength of selection

        if self.eq_z0:
            self.k = 3  # Number of parameters under the model
        else:
            self.k = 4

    def mean(self, test=False):
        if test:
            sd = np.diag(self.sd)[self.test_mask]
        else:
            sd = np.diag(self.sd)[self.training_mask]
        return self.z0 * np.exp(-self.alpha * sd) + (
                    1. - np.exp(-self.alpha * sd)) * self.z_eq

    def cov(self, test=False):
        if test:
            row_mask = self.test_mask
        else:
            row_mask = self.training_mask

        sd = self.sd[np.ix_(row_mask, self.training_mask)]
        pd = self.pd[np.ix_(row_mask, self.training_mask)]
        return np.exp(-self.alpha * pd) * (1. - np.exp(-2. * self.alpha * sd)) * (
                    self.sigma / (2. * self.alpha))

    def fit(self):
        # Fit the hyperparameters to the data

        init_params = np.random.uniform(size=self.k)

        bounds = [
            (None, None),
            (1e-12, None),
            (1e-12, None)
        ]

        if not self.eq_z0:
            bounds = [(None, None)] + bounds

        def objective(params):
            if self.eq_z0:
                self.z_eq, self.sigma, self.alpha = params
                self.z0 = self.z_eq
            else:
                self.z0, self.z_eq, self.sigma, self.alpha = params

            return self.nll()

        optim_res = optim.minimize(objective,
                                   init_params,
                                   bounds=bounds)

        return {
            'Optimization': {
                'Success': optim_res.success,
                'Message': optim_res.message
            },
            'Loglikelihood': -optim_res.fun,
            'DOF': self.k,
            'AIC': self.aic(optim_res.fun),
            'AIC.c': self.corrected_aic(optim_res.fun),
            'BIC': self.bic(optim_res.fun),
            'Parameters': dict(zip(['Z0', 'Zeq', 'sigma', 'alpha'][-self.k:],
                                   optim_res.x))
        }


class BM(GaussianProcessModel):

    def __init__(self, data, tree):
        super().__init__(data, tree)

        # The 4 parameters of the model:
        self.z0 = None  # The ancestral phenotype
        self.sigma = None  # The variance

        self.k = 2  # Number of parameters under the model

    def mean(self, test=False):
        return self.z0

    def cov(self, test=False):
        if test:
            row_mask = self.test_mask
        else:
            row_mask = self.training_mask

        sd = self.sd[np.ix_(row_mask, self.training_mask)]
        return sd * self.sigma

    def fit(self):
        # Fit the hyperparameters to the data

        init_params = np.random.uniform(size=self.k)

        bounds = [
            (None, None),
            (1e-12, None)
        ]

        def objective(params):
            self.z0, self.sigma = params
            return self.nll()

        optim_res = optim.minimize(objective,
                                   init_params,
                                   bounds=bounds)

        return {
            'Optimization': {
                'Success': optim_res.success,
                'Message': optim_res.message
            },
            'Loglikelihood': -optim_res.fun,
            'DOF': self.k,
            'AIC': self.aic(optim_res.fun),
            'AIC.c': self.corrected_aic(optim_res.fun),
            'BIC': self.bic(optim_res.fun),
            'Parameters': dict(zip(['Z0', 'sigma'],
                                   optim_res.x))
        }
