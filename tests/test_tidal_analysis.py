"""
Test that the function performs as expected with some synthetic signals
"""
from pathlib import Path
import sys
library_dir = Path(__file__).parent.parent.absolute()
print(library_dir)
sys.path.append(str(library_dir))
import math
import numpy as np
from quick_tidal_analysis.tidal_analysis import tidal_analysis


def test_sinlge_cosine():
    #Single arbitrary cosine accuracy requirements
    # amplitude to 1 part in 10000
    # phase to 1 degree
    periods = 12.43
    freqs = 24./periods
    amps= np.random.rand() +1.0 # minimum amplitude for detectability
    phas = np.pi * (np.random.rand() - 0.5)
    t = np.arange(0,60,0.001)
    x = amps*np.cos((2*np.pi*t*freqs)+phas) 
    __, const, __ =  tidal_analysis(x,t,T=[periods])

    assert math.isclose(const.amplitude.T1, amps, rel_tol=1e-4)
    assert math.isclose(const.phase.T1, phas, abs_tol=np.pi/180)



def test_out_phase_triple_cosine():
    # Triple out of phase radndom cosinses accuracy requirements:
    # amplitude to 1 part in 50
    # phase to 1 degree
    periods = np.array([12.43, 23.93, 16])
    freqs = 24/periods
    amps= np.random.rand(3) +1.0 # minimum amplitude for detectability
    phas = np.pi * (np.random.rand(3) - 0.5)
    t = np.arange(0,60,0.001)
    x = amps[0]*np.cos((2*np.pi*t*freqs[0])+phas[0]) + amps[1]*np.cos((2*np.pi*t*freqs[1])+phas[1])  + amps[2]*np.cos((2*np.pi*t*freqs[2])+phas[2])
    __, const, __ =  tidal_analysis(x,t,T=periods)
    for i in range(3):
        assert math.isclose(const.amplitude["T"+str(i+1)], amps[i], rel_tol=0.02)
        assert math.isclose(const.phase["T"+str(i+1)], phas[i], abs_tol=np.pi/180)

                                             
def test_noise():
    # Single noisy cosine wave accuracy requirements:
    # Noise is gaussian centered on 0 at half amplitude of wave
    # amp to 1 part in 10
    # phase to 2 degrees
    t = np.arange(0,20,0.01)
    phase =  np.pi * (np.random.rand() - 0.5)
    x = np.cos(t*2*np.pi + phase) + 0.5*np.random.randn((len(t)))
    __, const, __ =  tidal_analysis(x,t,T=[24])
    assert math.isclose(const.amplitude.T1, 1, rel_tol=0.05)
    assert math.isclose(const.phase.T1, phase, abs_tol=np.pi/90)

