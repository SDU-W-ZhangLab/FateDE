import sys
import logging
from time import time
import warnings
import matplotlib.pyplot as plt
#from plot import *
from numpy import polyfit
from scipy.stats.stats import pearsonr, spearmanr, kendalltau
import numpy as np
from scipy import optimize
from scipy import stats
from scipy.misc import derivative
import statsmodels.api as sm

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from tqdm.autonotebook import tqdm

import pandas as pd
import scipy as sp
from scipy import interpolate


def qvalue(pv, pi0=None):
    '''
    Estimates q-values from p-values

    This function is modified based on https://github.com/nfusi/qvalue

    Args
    ====
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.

    '''
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    m = float(len(pv))

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = sp.arange(0, 0.90, 0.01)
        counts = sp.array([(pv > i).sum() for i in sp.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = sp.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)

        if pi0 > 1:
            pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    p_ordered = sp.argsort(pv)
    pv = pv[p_ordered]
    qv = pi0 * m/len(pv) * pv
    qv[-1] = min(qv[-1], 1.0)

    for i in range(len(pv)-2, -1, -1):
        qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

    # reorder qvalues
    qv_temp = qv.copy()
    qv = sp.zeros_like(qv)
    qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv


def get_l_limits(X):
    Xsq = np.sum(np.square(X), 1) 
    R2 = -2. * np.dot(X, X.T) + (Xsq[:, None] + Xsq[None, :])
    R2 = np.clip(R2, 0, np.inf) 
    R_vals = np.unique(R2.flatten())
    R_vals = R_vals[R_vals > 1e-8]

    l_min = np.sqrt(R_vals.min()) / 2.
    l_max = np.sqrt(R_vals.max()) * 2.

    return l_min, l_max

## Kernels ##

def SE_kernel(X, l):
    Xsq = np.sum(np.square(X), 1)
    R2 = -2. * np.dot(X, X.T) + (Xsq[:, None] + Xsq[None, :])
    R2 = np.clip(R2, 1e-12, np.inf)
    return np.exp(-R2 / (2 * l ** 2))


def linear_kernel(X):
    K = np.dot(X, X.T)
    return K / K.max()


def cosine_kernel(X, p):
    ''' Periodic kernel as l -> oo in [Lloyd et al 2014]

    Easier interpretable composability with SE?
    '''
    Xsq = np.sum(np.square(X), 1)
    R2 = -2. * np.dot(X, X.T) + (Xsq[:, None] + Xsq[None, :])
    R2 = np.clip(R2, 1e-12, np.inf)
    return np.cos(2 * np.pi * np.sqrt(R2) / p)

def poly_kernel(X):
    
    K = np.power(np.dot(X,X.T), 2)
    return K

def sigmoid_kernel(X):
    

    K = np.tanh(0.01*np.dot(X, X.T))
    return K

def gower_scaling_factor(K):
    ''' Gower normalization factor for covariance matric K

    Based on https://github.com/PMBio/limix/blob/master/limix/utils/preprocess.py
    '''
    n = K.shape[0]
    P = np.eye(n) - np.ones((n, n)) / n
    KP = K - K.mean(0)[:, np.newaxis]
    trPKP = np.sum(P * KP)

    return trPKP / (n - 1)


def factor(K):
    S, U = np.linalg.eigh(K)
    # .clip removes negative eigenvalues
    return U, np.clip(S, 1e-8, None)


def get_UT1(U):
    return U.sum(0)


def get_UTy(U, y):
    return y.dot(U)


def mu_hat(delta, UTy, UT1, S, n, Yvar=None):
    ''' ML Estimate of bias mu, function of delta.
    '''
    if Yvar is None:
        Yvar = np.ones_like(S)

    UT1_scaled = UT1 / (S + delta * Yvar)
    sum_1 = UT1_scaled.dot(UTy)
    sum_2 = UT1_scaled.dot(UT1)

    return sum_1 / sum_2


def s2_t_hat(delta, UTy, S, n, Yvar=None):
    ''' ML Estimate of structured noise, function of delta
    '''
    if Yvar is None:
        Yvar = np.ones_like(S)

    UTy_scaled = UTy / (S + delta * Yvar)
    return UTy_scaled.dot(UTy) / n


def LL(delta, UTy, UT1, S, n, Yvar=None):
    ''' Log-likelihood of GP model as a function of delta.

    The parameter delta is the ratio s2_e / s2_t, where s2_e is the
    observation noise and s2_t is the noise explained by covariance
    in time or space.
    '''

    mu_h = mu_hat(delta, UTy, UT1, S, n, Yvar)

    if Yvar is None:
        Yvar = np.ones_like(S)

    sum_1 = (np.square(UTy - UT1 * mu_h) / (S + delta * Yvar)).sum()
    sum_2 = np.log(S + delta * Yvar).sum()

    with np.errstate(divide='ignore'):
        return -0.5 * (n * np.log(2 * np.pi) + n * np.log(sum_1 / n) + sum_2 + n)


def logdelta_prior_lpdf(log_delta):
    s2p = 100.
    return -np.log(np.sqrt(2 * np.pi * s2p)) - np.square(log_delta - 20.) / (2 * s2p)


def make_objective(UTy, UT1, S, n, Yvar=None):
    def LL_obj(log_delta):
        return -LL(np.exp(log_delta), UTy, UT1, S, n, Yvar)

    return LL_obj


def brent_max_LL(UTy, UT1, S, n):
    LL_obj = make_objective(UTy, UT1, S, n)
    o = optimize.minimize_scalar(LL_obj, bounds=[-10, 10], method='bounded', options={'maxiter': 32})
    max_ll = -o.fun
    max_delta = np.exp(o.x)
    max_mu_hat = mu_hat(max_delta, UTy, UT1, S, n)
    max_s2_t_hat = s2_t_hat(max_delta, UTy, S, n)

    return max_ll, max_delta, max_mu_hat, max_s2_t_hat


def lbfgsb_max_LL(UTy, UT1, S, n, Yvar=None):
    LL_obj = make_objective(UTy, UT1, S, n, Yvar)
    min_boundary = -10
    max_boundary = 20.
    x, f, d = optimize.fmin_l_bfgs_b(LL_obj, 0., approx_grad=True,
                                     bounds=[(min_boundary, max_boundary)],
                                     maxfun=64, factr=1e12, epsilon=1e-4)
    max_ll = -f
    max_delta = np.exp(x[0])

    boundary_ll = -LL_obj(max_boundary)
    if boundary_ll > max_ll:
        max_ll = boundary_ll
        max_delta = np.exp(max_boundary)

    boundary_ll = -LL_obj(min_boundary)
    if boundary_ll > max_ll:
        max_ll = boundary_ll
        max_delta = np.exp(min_boundary)

    max_mu_hat = mu_hat(max_delta, UTy, UT1, S, n, Yvar)
    max_s2_t_hat = s2_t_hat(max_delta, UTy, S, n, Yvar)

    s2_logdelta = 1. / (derivative(LL_obj, np.log(max_delta), n=2) ** 2)

    return max_ll, max_delta, max_mu_hat, max_s2_t_hat, s2_logdelta


def search_max_LL(UTy, UT1, S, n, num=32):
    ''' Search for delta which maximizes log likelihood.
    '''
    min_obj = np.inf
    max_log_delta = np.nan
    LL_obj = make_objective(UTy, UT1, S, n)
    for log_delta in np.linspace(start=-10, stop=20, num=num):
        cur_obj = LL_obj(log_delta)
        if cur_obj < min_obj:
            min_obj = cur_obj
            max_log_delta = log_delta

    max_delta = np.exp(max_log_delta)
    max_mu_hat = mu_hat(max_delta, UTy, UT1, S, n)
    max_s2_t_hat = s2_t_hat(max_delta, UTy, S, n)
    max_ll = -min_obj

    return max_ll, max_delta, max_mu_hat, max_s2_t_hat


def make_FSV(UTy, S, n, Gower):
    def FSV(log_delta):
        s2_t = s2_t_hat(np.exp(log_delta), UTy, S, n)
        s2_t_g = s2_t * Gower

        return s2_t_g / (s2_t_g + np.exp(log_delta) * s2_t)

    return FSV


def lengthscale_fits(exp_tab, U, UT1, S, Gower, num=64):
    ''' Fit GPs after pre-processing for particular lengthscale
    '''
    results = []
    n, G = exp_tab.shape
    for g in tqdm(range(G), leave=False):
        y = exp_tab.iloc[:, g]
        UTy = get_UTy(U, y)

        t0 = time()
        max_reg_ll, max_delta, max_mu_hat, max_s2_t_hat, s2_logdelta = lbfgsb_max_LL(UTy, UT1, S, n)
        max_ll = max_reg_ll
        t = time() - t0

        # Estimate standard error of Fraction Spatial Variance
        FSV = make_FSV(UTy, S, n, Gower)
        s2_FSV = derivative(FSV, np.log(max_delta), n=1) ** 2 * s2_logdelta

        results.append({
            'g': exp_tab.columns[g],  # gene name 
            'max_ll': max_ll,  # maximum likelihood
            'max_delta': max_delta,  # max delta
            'max_mu_hat': max_mu_hat,  # max mu
            'max_s2_t_hat': max_s2_t_hat,  # max sigma^2_t
            'time': t,
            'n': n,
            'FSV': FSV(np.log(max_delta)),  # FSV
            's2_FSV': s2_FSV,  # s^2_FSV
            's2_logdelta': s2_logdelta  # s^2_delta
        })

    return pd.DataFrame(results)


def null_fits(exp_tab):
    ''' Get maximum LL for null model
    '''
    results = []
    n, G = exp_tab.shape
    for g in range(G):
        y = exp_tab.iloc[:, g]
        max_mu_hat = 0.
        max_s2_e_hat = np.square(y).sum() / n  # mll estimate
        max_ll = -0.5 * (n * np.log(2 * np.pi) + n + n * np.log(max_s2_e_hat))

        results.append({
            'g': exp_tab.columns[g],
            'max_ll': max_ll,
            'max_delta': np.inf,
            'max_mu_hat': max_mu_hat,
            'max_s2_t_hat': 0.,
            'time': 0,
            'n': n
        })

    return pd.DataFrame(results)


def const_fits(exp_tab):
    ''' Get maximum LL for const model
    '''
    results = []
    n, G = exp_tab.shape
    for g in range(G):
        y = exp_tab.iloc[:, g]
        max_mu_hat = y.mean()
        max_s2_e_hat = y.var()
        sum1 = np.square(y - max_mu_hat).sum()
        max_ll = -0.5 * (n * np.log(max_s2_e_hat) + sum1 / max_s2_e_hat + n * np.log(2 * np.pi))

        results.append({
            'g': exp_tab.columns[g],
            'max_ll': max_ll,
            'max_delta': np.inf,
            'max_mu_hat': max_mu_hat,
            'max_s2_t_hat': 0.,
            'time': 0,
            'n': n
        })

    return pd.DataFrame(results)


def simulate_const_model(MLL_params, N):
    dfm = np.zeros((N, MLL_params.shape[0]))
    for i, params in enumerate(MLL_params.iterrows()):
        params = params[1]
        s2_e = params.max_s2_t_hat * params.max_delta
        dfm[:, i] = np.random.normal(params.max_mu_hat, s2_e, N)

    dfm = pd.DataFrame(dfm)
    return dfm


def get_mll_results(results, null_model='const'):
    null_lls = results.query('model == "{}"'.format(null_model))[['g', 'max_ll']]
    model_results = results.query('model != "{}"'.format(null_model))
    
    model_results = model_results[model_results.groupby(['g'])['max_ll'].transform(max) == model_results['max_ll']]
    mll_results = model_results.merge(null_lls, on='g', suffixes=('', '_null'))
    mll_results['LLR'] = mll_results['max_ll'] - mll_results['max_ll_null']

    return mll_results


def dyn_de(X, exp_tab, kernel_space=None):
    if kernel_space == None:
        kernel_space = {
            'SE': [5., 25., 50.]
        }

    results = []

    if 'null' in kernel_space:
        result = null_fits(exp_tab)
        result['l'] = np.nan
        result['M'] = 1
        result['model'] = 'null'
        results.append(result)

    if 'const' in kernel_space:
        result = const_fits(exp_tab)  
        result['l'] = np.nan
        result['M'] = 2
        result['model'] = 'const'
        results.append(result)

    logging.info('Pre-calculating USU^T = K\'s ...')
    US_mats = []
    t0 = time()
    if 'linear' in kernel_space:
        K = linear_kernel(X)
        U, S = factor(K)
        gower = gower_scaling_factor(K)
        UT1 = get_UT1(U)
        US_mats.append({
            'model': 'linear',
            'M': 3,
            'l': np.nan,
            'U': U,
            'S': S,
            'UT1': UT1,
            'Gower': gower
        })

    if 'poly' in kernel_space:
        K = poly_kernel(X)
        U, S = factor(K)
        gower = gower_scaling_factor(K)
        UT1 = get_UT1(U)
        US_mats.append({
            'model': 'poly',
            'M': 3,
            'l': np.nan,
            'U': U,
            'S': S,
            'UT1': UT1,
            'Gower': gower
        })
    if 'SE' in kernel_space:
        for lengthscale in kernel_space['SE']:
            K = SE_kernel(X, lengthscale)  
            U, S = factor(K)  
            gower = gower_scaling_factor(K)  # the Gower factor for covariance matrix
            UT1 = get_UT1(U)  
            US_mats.append({
                'model': 'SE',
                'M': 4,
                'l': lengthscale,
                'U': U,
                'S': S,
                'UT1': UT1,
                'Gower': gower
            })

    if 'PER' in kernel_space:
        for period in kernel_space['PER']:
            K = cosine_kernel(X, period)
            U, S = factor(K)
            gower = gower_scaling_factor(K)
            UT1 = get_UT1(U)
            US_mats.append({
                'model': 'PER',
                'M': 4,
                'l': period,
                'U': U,
                'S': S,
                'UT1': UT1,
                'Gower': gower
            })

    t = time() - t0
    logging.info('Done: {0:.2}s'.format(t))

    logging.info('Fitting gene models')
    n_models = len(US_mats)
    for i, cov in enumerate(tqdm(US_mats, desc='Models: ')):
        result = lengthscale_fits(exp_tab, cov['U'], cov['UT1'], cov['S'], cov['Gower'])  
        result['l'] = cov['l']
        result['M'] = cov['M']
        result['model'] = cov['model']
        results.append(result)

    n_genes = exp_tab.shape[1]
    logging.info('Finished fitting {} models to {} genes'.format(n_models, n_genes))
    results = pd.concat(results, sort=True).reset_index(drop=True)
    results['BIC'] = -2 * results['max_ll'] + results['M'] * np.log(results['n'])

    return results

def run_per(X, exp_tab, kernel_space=None):
    ''' Perform SpatialDE test

        X : matrix of spatial coordinates times observations
        exp_tab : Expression table, assumed appropriatealy normalised.

        The grid of covariance matrices to search over for the alternative
        model can be specifiec using the kernel_space paramter.
        '''
    P_min, P_max = 0., 0.
    if kernel_space == None:
        P_min, P_max = get_l_limits(X)

        kernel_space = {
            'PER': np.logspace(np.log10(P_max), np.log10(P_max*4.), 2),
            'const': 0
        }
    logging.info('Performing PER test')
    results = dyn_de(X, exp_tab, kernel_space)
    mll_results = get_mll_results(results)  
    # Perform significance test
    mll_results['pval'] = 1 - stats.chi2.cdf(mll_results['LLR'], df=1)  
    mll_results['qval'] = qvalue(mll_results['pval'])  

    return mll_results, P_min, P_max

def logNorm(df):
    dfN = df
    for i in df.columns:
        if i != "Pseudotime"  and i != "State":
            dfN[i] = np.log(df[i])
    return dfN

def G_result(results_gene, fate_name, P_max):
    G1_1 = []
    for index, row in results_gene.iterrows():
        pair_tmp = {
            "Branch": fate_name,
            "gene": row['g'],
            "l": row['l'] / P_max,
            "pval": row['pval'],
            "FSV": row['FSV']
        }
        G1_1.append(pair_tmp)
    G1_1 = pd.DataFrame(G1_1).sort_values(by=['l', 'FSV'], ascending=[True, False])
    return G1_1



def find_genes(df1, df2, G1, G2, fate1_name, fate2_name,logN=True):
    """Find genes that increase first and then decrease in another branch
        df1,df2：DataFrame,
            Expression matrix of two branches, behavioral cells, listed as genes, with Pseudotime column
        G1, G2: list
        fate1_name, fate2_name: str, fate name
        logN：if log
    """
    row, col = df1.shape
    if logN:
        print("Log Normalization ............ ")
        df1 = logNorm(df1)
        df2 = logNorm(df2)

    df1 = df1.sort_values(by='Pseudotime')
    df2 = df2.sort_values(by='Pseudotime')

    # Gaussian process
    X = df2[['Pseudotime']]  
    Y = df2.loc[:, G1] 

    results1_gene, _, P1_max = run_per(pd.DataFrame(X), Y)
    G1_1 = G_result(results1_gene, fate1_name, P1_max)

    # Gaussian process
    X = df1[['Pseudotime']]  
    Y = df1.loc[:, G2]  

    results2_gene, _, P2_max = run_per(pd.DataFrame(X), Y)
    G2_1 = G_result(results2_gene, fate2_name, P2_max)

    return G1_1, G2_1

# def fitting(df, fit_degree):
#     df_fit = df
#     X_array = np.array(df_fit.loc[:, "Pseudotime"].values.tolist())
#     for gene in df_fit.columns:
#         if gene != "Pseudotime"  and gene != "State": 
#             Y_array = np.array(df_fit.loc[:, gene].values.tolist())
#             coeff = polyfit(X_array, Y_array, fit_degree)
#             _y = np.polyval(coeff, X_array)
#             df_fit[gene] = _y
#     return df_fit
# 
# def cal_pair(df, G1_1, G2_1):
#     fate1_tmp = df.loc[:, G1_1 + G2_1 + ["Pseudotime"]]
#     fate1_tmp = fate1_tmp.diff().dropna()
#     fate1_diff = pd.DataFrame()
#     for i in range(fate1_tmp.shape[1] - 1):
#         gene = fate1_tmp.columns[i]
#         # fate1_diff[gene] = fate1_tmp.iloc[:,i]/fate1_tmp.iloc[:, fate1_tmp.shape[1]-1]
#         fate1_diff[gene] = fate1_tmp.iloc[:, i]  
#     return fate1_diff

# def find_pair(df1, df2, G1_1, G2_1, State=[1,2,3], fit_degree=0,corr="pearsonr"):
#     """ find pair
#         df: exp matrix
#         State:list
#         fit_degree: 
#         corr:  pearsonr, spearmanr, kendalltau
#     """
#     if fit_degree > 0:
#         print("Fitting .......")
#         df1 = fitting(df1, fit_degree)
#         df2 = fitting(df2, fit_degree)
# 
#     df1 = df1.sort_values(by='Pseudotime')
#     df2 = df2.sort_values(by='Pseudotime')
#     if df1.shape[0] != df2.shape[0]:
#         print("uniform......")
#         df1_state = df1[df1.State == State[1]]
#         df2_state = df2[df2.State == State[2]]
#         min_i = min(df1_state.shape[0], df2_state.shape[0])
# 
#         df1_index = np.linspace(0, df1_state.shape[0] - 1, min_i, dtype=int)
#         df2_index = np.linspace(0, df2_state.shape[0] - 1, min_i, dtype=int)  # get index
#         df1_state = df1_state.iloc[df1_index]
#         df2_state = df2_state.iloc[df2_index]  
# 
#         df1 = pd.concat([df1[df1.State == State[0]], df1_state])
#         df2 = pd.concat([df2[df2.State == State[0]], df2_state])
# 
#     # Create derivative of Fate 1
#     fate1_diff = cal_pair(df1, G1_1, G2_1)
#     # print(fate1_diff)
#     # Create a pair for Fate 1
#     fate1_pair_result = pd.DataFrame()
#     for dx2 in G2_1:
#         for dx1 in G1_1:
#             dx2_dx1 = fate1_diff.loc[:, dx2] / fate1_diff.loc[:, dx1].add(0.00001)  
#             fate1_pair_result[dx1 + "_" + dx2] = dx2_dx1  # Pay attention to naming
# 
#     # Create derivative of Fate 2
#     fate2_diff = cal_pair(df2, G1_1, G2_1)
#     # Create a pair for Fate 2
#     fate2_pair_result = pd.DataFrame()
#     for dx2 in G2_1:
#         for dx1 in G1_1:
#             dx1_dx2 = fate2_diff.loc[:, dx1] / fate2_diff.loc[:, dx2].add(0.00001)
#             fate2_pair_result[dx1 + "_" + dx2] = dx1_dx2  # Pay attention to naming
# 
#     # Ensure that the number of rows taken is the same
#     min_cell_num = min(fate1_pair_result.shape[0], fate2_pair_result.shape[0])
#     # fate1_pair_result = fate1_pair_result.iloc[:min_cell_num, :]
#     # fate2_pair_result = fate2_pair_result.iloc[:min_cell_num, :]
#     fate1_pair_result.index = range(min_cell_num)
#     fate2_pair_result.index = range(min_cell_num)
# 
#     
#     fate_pair_result = []
# 
#     if corr == "pearsonr":
#         for pair in fate1_pair_result.columns:
#             r, p = pearsonr(fate1_pair_result.loc[:, pair].tolist(), fate2_pair_result.loc[:, pair].tolist())
#             # est = sm.OLS(fate1_pair_result.loc[:, pair], fate2_pair_result.loc[:, pair])
#             # est2 = est.fit()
#             fate_pair_result.append({"pair": pair, "r": r, "pvalues":p})
#     elif corr == "spearmanr":
#         for pair in fate1_pair_result.columns:
#             r, p = spearmanr(fate1_pair_result.loc[:, pair].tolist(), fate2_pair_result.loc[:, pair].tolist())
# 
#             fate_pair_result.append({"pair": pair, "r": r, "pvalues": p})
#     elif corr == "kendalltau":
#         for pair in fate1_pair_result.columns:
#             r, p = kendalltau(fate1_pair_result.loc[:, pair].tolist(), fate2_pair_result.loc[:, pair].tolist())
#             fate_pair_result.append({"pair": pair, "r": r, "pvalues": p})
#     fate_pair_result = pd.DataFrame(fate_pair_result).sort_values(by="r", ascending=False)
# 
#     return fate_pair_result, df1, df2,fate1_pair_result,fate2_pair_result
