__doc__ = """\
a module holding simple utilities for reading bayesline chains. The conventions assumed herein are taken from the BayesWave source code:
    spline file structure: BayesWaveIO.c line 1737
    spline model:          BayesLine.c lines 1105-1139
    line file structure:   BayesWaveIO.c line 1773
    line model:            BayesWaveIO.c lines 20-66
"""
__author__ = "reed.essick@ligo.org"

#-------------------------------------------------

import numpy as np
from scipy.interpolate import UnivariateSpline

#-------------------------------------------------

DEFAULT_BURNIN = 0
DEFAULT_DOWNSAMPLE = 1

SPLINE_DTYPE = [('frequency', 'float'), ('lnPSD', 'float')]
LINES_DTYPE = [('frequency', 'float'), ('amplitude', 'float'), ('quality', 'float')]

#-------------------------------------------------

def splinepath2samples(path, burnin=DEFAULT_BURNIN, downsample=DEFAULT_DOWNSAMPLE):
    """
    reads in samples from file, discards the first "burnin" samples and then retains one out of every "downsample" remaining samples
    """
    samples = []
    with open(path, 'r') as file_obj:
        for ind, line in enumerate(file_obj):
            if (ind>=burnin) and (ind%downsample==0):
                fields = [float(_) for _ in line.strip().split()[1:]]
                samples.append(np.array(zip(fields[::2], fields[1::2]), dtype=SPLINE_DTYPE)) ### stored as (frequency, ln(PSD)) pairs

    return samples

def linespath2samples(path, burnin=DEFAULT_BURNIN, downsample=DEFAULT_DOWNSAMPLE):
    """
    reads in samples from file, discards the first "burnin" samples and then retains one out of every "downsample" remaining samples
    """
    samples = []
    with open(path, 'r') as file_obj:
        for ind, line in enumerate(file_obj):
            if (ind>=burnin) and (ind%downsample==0):
                fields = [float(_) for _ in line.strip().split()[1:]]
                samples.append(np.array(zip(fields[::3], fields[1::3], fields[2::3]), dtype=LINES_DTYPE)) ### stored as (frequency, ln(PSD)) pairs

    return samples

def samples2psd(freqs, spline_samples, lines_samples, verbose=False):
    """
    returns psd estimates based on the samples provided
    """
    assert len(spline_samples)==len(lines_samples), 'we must have the same number of spline_samples and line_samples'
    return splinesamples2psd(freqs, spline_samples, verbose=verbose) + linessamples2psd(freqs, lines_samples, verbose=verbose)

def splinesamples2psd(freqs, samples, verbose=False):
    """
    return the spline's contribution to the PSD estimate at freqs for each sample
    """
    nfreq = len(freqs)
    nsamp = len(samples)
    psds = np.empty((nsamp, nfreq), dtype=float)
    for ind, sample in enumerate(samples):
        if verbose:
            print('spline sample %d'%ind)
        spline = UnivariateSpline(sample['frequency'], sample['lnPSD'], k=3, s=0) ### k=3 --> cubic, s=0 --> interpolate through data
        psds[ind,:] = np.exp(spline(freqs))
    
    return psds

def linessamples2psd(freqs, samples, verbose=False):
    """
    return the lines' contribution to the PSD estimate at freqs for each sample
    """
    nfreq = len(freqs)
    nsamp = len(samples)
    psds = np.zeros((nsamp, nfreq), dtype=float)
    for ind, sample in enumerate(samples):
        if verbose:
            print('line sample %d'%ind)
        for f, a, q in sample:
            psds[ind,:] += line(freqs, f, a, q)

    return psds

def line(freqs, f, a, q):
    """
    the line model Bayesline implements
    """
    psd = np.zeros_like(freqs, dtype=float)

    ### figure out which frequencies actually contribute
    spread = max(50, q*1e-2)
    df = f/spread
    truth = (f-df <= freqs)*(freqs <= f+df)
    freqs = freqs[truth]

    ### compute contribution
    z = np.exp(-(np.abs(freqs-f)-df)/df)
    z[z>1] = 1

    f2 = f**2
    freqs2 = freqs**2

    psd[truth] = z*(a*f2/q**2)/(freqs2*(freqs2-f2)**2)
    return psd
