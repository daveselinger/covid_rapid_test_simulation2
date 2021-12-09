
from scipy.stats import norm

def gaussianRandom(mu,var=1.0):
    return norm.rvs(mu,var)

