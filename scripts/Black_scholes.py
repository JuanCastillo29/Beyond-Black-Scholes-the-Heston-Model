# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 11:07:41 2025

@author: juanc
"""

import numpy as np
from scipy.stats import norm

def Black_Scholes_pricing(T, K, S, r_free, vol, option_type='call'):
    d1 = (np.log(S/K) + (r_free + 0.5 * vol**2) * T) / (vol * np.sqrt(T))
    d2 = d1 - vol * np.sqrt(T)
    if option_type == 'call':
        return S * norm.cdf(d1) - K * np.exp(-r_free * T) * norm.cdf(d2)
    else:
        return K * np.exp(-r_free * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
    
def Black_Scholes_Delta(T, K, S, r_free, vol, option_type='call'):
    d1 = (np.log(S/K) + (r_free + 0.5 * vol**2) * T) / (vol * np.sqrt(T))
    
    if option_type == 'call':
        return norm.cdf(d1)
    else:
        return norm.cdf(d1) - 1
    
def Black_Scholes_Gamma(T, K, S, r_free, vol):
    d1_mod = (np.log(S/K) + (r_free + 0.5 * vol**2) * T) 
    d1_2 = d1_mod*d1_mod / (vol**2 * T)
    return np.exp(- 0.5 * d1_2)/np.sqrt(2 * np.pi * T) / (S * vol)

def Black_Scholes_Vega(T, K, S,Gamma, vol):
    return Gamma * S * S * vol * T

    
def Black_Scholes_Rho(T, K, S, r_free, vol, option_type='call'):
    d1 = (np.log(S/K) + (r_free + 0.5 * vol**2) * T) / (vol * np.sqrt(T))
    d2 = d1 - vol * np.sqrt(T)
    rho = K * T * np.exp(-r_free * T)
    if option_type == 'call':
        return rho * norm.cdf(d2)
    else:
        return rho * (norm.cdf(d2) - 1.0)
    
def Black_Scholes_Theta(T, K, S, r_free, vol,Gamma, Rho):
    return - S * S * vol * vol * Gamma/2.0 - Rho * r_free