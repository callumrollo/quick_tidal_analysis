# Quick tidal analysis

A python package to extract harmonic components from a time series. 

Try it out from the comfort of your browser: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/callumrollo/quick_tidal_analysis/main?filepath=Explore.ipynb)

May take a minute to load

### Under development

This is not a production ready script, It is still a work in progress with no tests yet.

Once more mature I might put it on PyPI and conda-forge

### Usage

To extract periodic sinusoids at 12 and 24 hour periods from a variable `x` sampled at times `t` in yeardays and plot the results

`df, df_consts, df_fine = tidal_analysis(x,t,T=[12,24],  plot_results=True)`

##### Returns:
`df` pandas dataframe indexed to input times with original x, detrended x, tidal sinusoids and residual x\
`df_const` pandas dataframe of tidal constituents' amplitude, phase, period and mean\
`df_fine` finely sampled (0.001 days) time series of tidal signals\
    

Check out the notebook [Explore.ipynb](https://github.com/callumrollo/quick_tidal_analysis/blob/main/Explore.ipynb) for more examples
