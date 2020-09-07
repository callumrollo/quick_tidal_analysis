"""
Quick tidal analysis script

Created by Callum Rollo
2020-09-03
"""

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import detrend
from scipy.optimize import leastsq

def tidal_analysis(x, t, T=[12.4206], remove_trend=True, plot_results=False, freq_seek=False):
    """
    This function performs a harmonic analysis of timeseries x for a single
    or multiple tidal constituents including removal of trend and time-mean.
    :param x: variable (any units)
    :param t: time (yearday)
    :param T: tidal periods (hours) [T1,T2,T3,...] default: 12.4206 (M2)
    :remove_trend: switch to remove trend as well as time-mean. Default=True
    :freq_seek: if True, cost function will seek in frequency space as well as amp and phase. Default=False
    :plot_results: create a plot of the signal and fitted sine curves. Default=False
    :returns:
    
    df: pandas dataframe indexed to input times with original x, detrended x, tidal sinusoids and residual x
    df_const: pandas dataframe of tidal constituents' amplitude, phase, period and mean
    df_fine: finely sampled (0.001 days) time series of tidal signals
    
    Original matlab function by Rob Hall (Aug 2014)
    Adapted to Python by Callum Rollo (Sep 2020)
    """  
    # Remove infinite and nan values from input
    t = t[np.isfinite(x)]
    x = x[np.isfinite(x)]

    # Create dataframe for variables
    df = pd.DataFrame({"time_yday": t, "x": x})
    df = df.set_index("time_yday")

    # Convert tidal periods to days
    T = np.array(T)
    T_days = T/24

    # subtract mean and detrend
    x_mean = np.nanmean(x)

    if remove_trend:
        x_proc = detrend(x - x_mean)
    else:
        x_proc = x - x_mean

    df['x_detrend'] = x_proc
    df['T_summed'] = 0
    
    # Finescale time for later upsampling of tidal signals
    t_fine = np.arange(t[0], t[-1]+0.001, 0.001)
    df_fine = pd.DataFrame({"time_yday":t_fine})
    df_fine['T_summed'] = 0

    # First guess values for the sinusoid amplitude based on signal variability
    guess_amp = np.std(x_proc)#/len(T_days)

    # Empty lists for tidal constituent info
    const_ids, amps, periods, phases, means = [], [], [], [], []
    
    # Loop through the supplied tidal periods T
    for i, time_period in enumerate(T_days):
        # Start from the frequency provided by user
        guess_freq = 2*np.pi/time_period

        # Define the function to optimize
        if freq_seek:
            # Optimises a cosine function by adjusting amp, frequency, phase and offset
            optimize_func = lambda x: x[0]*np.cos((x[1]*t)+x[2]) + x[3] - x_proc
            est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, 0, 0])[0]
            # Can return phases with mag greater than pi for some reason
            while np.abs(est_phase)>np.pi:
                est_phase = (np.abs(est_phase) - 2*np.pi) * est_phase/np.abs(est_phase)
        else:
            # Optimises a cosine function by adjusting amp, phase and offset
            optimize_func = lambda x: x[0]*np.cos(guess_freq*t+x[1]) + x[2] - x_proc
            est_amp,  est_phase, est_mean = leastsq(optimize_func, [guess_amp, 0, 0])[0]  
            est_freq = guess_freq
        # Append the values defining the cosine curve to the lists
        amps.append(np.abs(est_amp))
        periods.append(24*2*np.pi/est_freq)
        phases.append(est_phase)
        means.append(est_mean)
        const_ids.append("T"+str(i+1))
        
        # Recreare the fitted cosine curve using the optimized parameters
        data_fit = est_amp * np.cos(est_freq*t+est_phase) + est_mean
        df[f"x_T{i+1}"] = data_fit
        df['T_summed'] = df['T_summed'] + data_fit
        
        # upsample the sine waves for plotting
        df_fine[f"x_T{i+1}"] = est_amp * np.cos(est_freq*t_fine+est_phase) + est_mean   
        df_fine['T_summed'] = df_fine['T_summed'] + est_amp * np.cos(est_freq*t_fine+est_phase) + est_mean   
        
    # Calculate the residual by subtracting the combined tidal signal
    df['x_residual'] = df.x - df.T_summed
    
    # Create a data frame for tidal constituents
    df_consts = pd.DataFrame({"const_id": const_ids, "period_hours": T, "amplitude":amps, "fitted_period_hours":periods, "phase":phases, "mean":means})
    df_consts = df_consts.set_index('const_id')   
    
    # Plot results of calculations
    if plot_results:
        num_plots = 2+len(T)
        fig, ax = plt.subplots(math.ceil(num_plots/2),2, figsize=(20,14), sharex=True)
        ax = ax.ravel()
        ax[0].plot(t, x, label="Original data")
        ax[0].plot(t,  df['x_residual'], label="Residual data")
        for const in range(len(T)):
            ax[const+1].plot(t,  df['x_detrend'], label='detrended data')
            ax[const+1].plot(t_fine, df_fine[f"x_T{const+1}"], label=f"sine period {np.round(df_consts['fitted_period_hours'][const],2)} hours")
        ax[const+2].plot(t,  df['x_detrend'], label='detrended data')
        ax[const+2].plot(t_fine,  df_fine['T_summed'], label='Summed constituents')
        for axis in ax[:const+3]:
            axis.legend()
        ax[-2].set(xlabel='Time (yeardays)')
        ax[-1].set(xlabel='Time (yeardays)')
 
    return df, df_consts, df_fine
