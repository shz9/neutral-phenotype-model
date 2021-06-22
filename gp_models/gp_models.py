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

    def __init__(self, data, tree, loss='nll'):

        self.tree = tree
        self.data = data
        self.loss = loss

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

    def pairwise_divergence_loss(self):
        """
        NOTE: This loss function is non-convex and we do not recommend
        fitting models using it. Included only for experimental purposes.
        """

        obs_phenotypes = self.data[np.array(self.tips)[self.training_mask]]
        exp_mean = self.mean()
        exp_covariance = self.cov()

        # Pairwise divergence loss:
        obs_div = []  # Observed divergence
        exp_div = []  # Expected divergence

        for i in range(len(obs_phenotypes)):
            for j in range(i + 1, len(obs_phenotypes)):
                obs_div.append((obs_phenotypes[i] - obs_phenotypes[j])**2)
                exp_div.append(
                        exp_covariance[i, i] + exp_covariance[j, j] - 2.*exp_covariance[i, j] +
                        (exp_mean[i] - exp_mean[j])**2
                )

        obs_div = np.array(obs_div)
        exp_div = np.array(exp_div)

        # Mean squared loss between observed divergence
        # and expected divergence:
        return (np.abs(obs_div - exp_div)).mean()

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
            print(e)
            return np.inf

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

    def __init__(self, data, tree, equilibrium_z0=False, fixed_params=None):

        super().__init__(data, tree)

        self.fixed_params = fixed_params or {}

        self.eq_z0 = equilibrium_z0

        # The 5 parameters of the model:
        self.z0 = None  # The ancestral phenotype
        self.z_eq = None  # The equilibrium value of the phenotype
        self.psi = None  # Psi, captures the relationship between ancestral and equilibrium sequence
        self.sigma_eq = None  # The equilibrium variance
        self.u = None  # The scaling factor u, or the rate in our context

        self.k = 5 - len(self.fixed_params)  # Number of parameters under the model
        if self.eq_z0:
            if 'Psi' not in self.fixed_params:
                self.k -= 1
            if 'Zeq' not in self.fixed_params:
                self.k -= 1

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

        data_mean = self.data[np.array(self.tips)[self.training_mask]].mean()

        if self.eq_z0:
            init_params = np.array([
                data_mean,
                np.random.uniform(),
                np.random.uniform()
            ])

            bounds = [
                (None, None),
                (1e-12, 1e12),
                (1., 2.)
            ]

        else:

            init_params = np.array([
                data_mean,
                data_mean,
                np.random.uniform(),
                np.random.uniform(),
                np.random.uniform()
            ])

            bounds = [
                (None, None),
                (None, None),
                (1e-12, 1e12),
                (1e-12, 1e12),
                (1., 2.)
            ]

        def objective(params):

            if self.eq_z0:
                self.z0, self.sigma_eq, self.u = params
                self.z_eq = self.z0
                self.psi = 2.*self.sigma_eq
            else:
                self.z0, self.z_eq, self.psi, self.sigma_eq, self.u = params

            for k, v in self.fixed_params.items():
                if k == 'Z0':
                    if self.eq_z0:
                        self.z_eq = self.z0 = v
                    else:
                        self.z0 = v
                elif k == 'Zeq':
                    if self.eq_z0:
                        self.z_eq = self.z0 = v
                    else:
                        self.z_eq = v
                elif k == 'Psi':
                    self.psi = v
                elif k == 'sigma_eq':
                    if self.eq_z0:
                        self.sigma_eq = v
                        self.psi = 2.*v
                    else:
                        self.sigma_eq = v
                elif k == 'u':
                    self.u = v

            if self.loss == 'nll':
                return self.nll()
            else:
                return self.pairwise_divergence_loss()

        optim_res = optim.minimize(objective,
                                   init_params,
                                   method='L-BFGS-B',
                                   bounds=bounds,
                                   options={'maxiter': 100})

        if self.eq_z0:
            inf_params = dict(zip(['Z0', 'sigma_eq', 'u'],
                                  optim_res.x))
            inf_params['Zeq'] = inf_params['Z0']
            inf_params['Psi'] = 2.*inf_params['sigma_eq']
        else:
            inf_params = dict(zip(['Z0', 'Zeq', 'Psi', 'sigma_eq', 'u'],
                                  optim_res.x))

        inf_params.update(self.fixed_params)

        if self.loss == 'nll':
            nll = optim_res.fun
            div = self.pairwise_divergence_loss()
        else:
            div = optim_res.fun
            nll = self.nll()

        return {
            'Optimization': {
                'Success': optim_res.success,
                'Message': optim_res.message
            },
            'Loglikelihood': -nll,
            'Pairwise divergence loss': div,
            'DOF': self.k,
            'AIC': self.aic(nll),
            'AIC.c': self.corrected_aic(nll),
            'BIC': self.bic(nll),
            'Parameters': inf_params
        }


class OU(GaussianProcessModel):

    def __init__(self, data, tree, equilibrium_z0=False, init_BM=True):

        super().__init__(data, tree)

        self.bm_params = None

        # If the init_BM flag is set to true,
        # we initialize the fitting procedure with
        # the maximum-likelihood estimates of the BM model parameters.
        if init_BM:
            bm = BM(data, tree)
            self.bm_params = bm.fit()['Parameters']

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

        if self.bm_params is not None:
            # If the BM model parameters are provided,
            # initialize corresponding parameters with
            # those estimates.
            if self.eq_z0:
                init_params = np.array([self.bm_params['Z0'],
                                        self.bm_params['sigma'], .1])
            else:
                init_params = np.array([self.bm_params['Z0'], self.bm_params['Z0'],
                                        self.bm_params['sigma'], .1])
        else:
            # Otherwise, initialize randomly.
            data_mean = self.data[np.array(self.tips)[self.training_mask]].mean()
            if self.eq_z0:
                init_params = [data_mean, np.random.uniform(),
                               np.random.uniform()]
            else:

                init_params = [data_mean, data_mean,
                               np.random.uniform(), np.random.uniform()]

            init_params = np.array(init_params)

        bounds = [
            (None, None),
            (1e-12, 1e12),
            (1e-12, 1e12)
        ]

        if not self.eq_z0:
            bounds = [(None, None)] + bounds

        def objective(params):
            if self.eq_z0:
                self.z_eq, self.sigma, self.alpha = params
                self.z0 = self.z_eq
            else:
                self.z0, self.z_eq, self.sigma, self.alpha = params

            if self.loss == 'nll':
                return self.nll()
            else:
                return self.pairwise_divergence_loss()

        optim_res = optim.minimize(objective,
                                   init_params,
                                   method='L-BFGS-B',
                                   bounds=bounds,
                                   options={'maxiter': 100})

        if self.loss == 'nll':
            nll = optim_res.fun
            div = self.pairwise_divergence_loss()
        else:
            div = optim_res.fun
            nll = self.nll()

        return {
            'Optimization': {
                'Success': optim_res.success,
                'Message': optim_res.message
            },
            'Loglikelihood': -nll,
            'Pairwise divergence loss': div,
            'DOF': self.k,
            'AIC': self.aic(nll),
            'AIC.c': self.corrected_aic(nll),
            'BIC': self.bic(nll),
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
        if test:
            return np.repeat(self.z0, len(self.test_mask))
        else:
            return np.repeat(self.z0, len(self.training_mask))

    def cov(self, test=False):
        if test:
            row_mask = self.test_mask
        else:
            row_mask = self.training_mask

        sd = self.sd[np.ix_(row_mask, self.training_mask)]
        return sd * self.sigma

    def fit(self):
        # Fit the hyperparameters to the data

        data_points = self.data[np.array(self.tips)[self.training_mask]]
        init_params = np.array([data_points.mean(),
                                np.random.uniform()])

        bounds = [
            (None, None),
            (1e-12, 1e12)
        ]

        def objective(params):
            self.z0, self.sigma = params

            if self.loss == 'nll':
                return self.nll()
            else:
                return self.pairwise_divergence_loss()

        optim_res = optim.minimize(objective,
                                   init_params,
                                   method='L-BFGS-B',
                                   bounds=bounds,
                                   options={'maxiter': 100})

        if self.loss == 'nll':
            nll = optim_res.fun
            div = self.pairwise_divergence_loss()
        else:
            div = optim_res.fun
            nll = self.nll()

        return {
            'Optimization': {
                'Success': optim_res.success,
                'Message': optim_res.message
            },
            'Loglikelihood': -nll,
            'Pairwise divergence loss': div,
            'DOF': self.k,
            'AIC': self.aic(nll),
            'AIC.c': self.corrected_aic(nll),
            'BIC': self.bic(nll),
            'Parameters': dict(zip(['Z0', 'sigma'],
                                   optim_res.x))
        }
