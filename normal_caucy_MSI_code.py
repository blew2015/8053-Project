import numpy as np
import scipy
from scipy import stats
import pandas as pd
import pymc as pm
#import theano.tensor as tt
from scipy import optimize
# matplotlib for plotting
import arviz as az
import time, sys, os



ind = int(sys.argv[1])

RANDOM_SEED = 8053
np.random.seed(RANDOM_SEED)

alpha = [0, 0, -10, -10, 10, 10]
beta = [1, 100, 1, 100, 1, 100]


theta_posteriorMean = []

for index in range(1000):
    data =[]
    for i in range(10):
        data.append(np.random.normal(0, 1, 1)[0])

    N_SAMPLES = 5000

    with pm.Model() as morph_model:
        theta = pm.Cauchy('theta', alpha=alpha[ind], beta = beta[ind])
    
        observed = pm.Normal('obs', mu=theta, sigma=1, observed=data)
        start=pm.find_MAP()


    
        # Using Metropolis Hastings Sampling
        step = pm.Metropolis()
    
        # Sample from the posterior using the sampling method
        lya_trace = pm.sample(N_SAMPLES, step=step, return_inferencedata=True)
        
    pos = az.summary(lya_trace,var_names=['theta'])
    theta_posteriorMean.append(pos['mean'][0])

output = {'name': theta_posteriorMean}

df = pd.DataFrame(output)
    
name = 'output10Data_' + str(alpha[ind]) + '_' + str(beta[ind]) + '.txt'

df.to_csv(name, sep='\t', index=False)
