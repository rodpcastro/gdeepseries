import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

eps = np.finfo(np.float64).eps

pltdir = 'error/'
os.makedirs(pltdir, exist_ok=True)

def plot_error(X, Y, ZI, ZS, title='', fname=''):
    abs_err_Z = np.abs(ZI - ZS)
    abs_err_Z[abs_err_Z < eps] = eps
    min_abs_err_Z = eps
    max_abs_err_Z = abs_err_Z.max()

    fig, ax = plt.subplots(figsize=(7,6))
    
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    
    lev_exp = np.arange(np.floor(np.log10(min_abs_err_Z)),
                        np.ceil(np.log10(max_abs_err_Z))+1)
    levs = np.power(10, lev_exp)
    
    cs = ax.contourf(X, Y, abs_err_Z, levs, locator=ticker.LogLocator())
    
    cbar = fig.colorbar(cs)
    
    fig.savefig(fname+'.svg', bbox_inches='tight')


if __name__ == '__main__':
    x = np.load('scipy/x.npy')
    y = np.load('scipy/y.npy')

    X, Y = np.meshgrid(x, y, indexing='ij')

    scp_time = np.load('scipy/time.npy')
    scp_f = np.load('scipy/f.npy')
    scp_fx = np.load('scipy/fx.npy')
    scp_fxx = np.load('scipy/fxx.npy')

    gds_time = np.load('gds/time.npy')
    gds_f = np.load('gds/f.npy').T
    gds_fx = np.load('gds/fx.npy').T
    gds_fxx = np.load('gds/fxx.npy').T

    plot_error(X, Y, scp_f, gds_f, title=r'$F$ absolute error', fname=pltdir+'errorf')
    plot_error(X, Y, scp_fx, gds_fx, title=r'$F_x$ absolute error', fname=pltdir+'errorfx')
    plot_error(X, Y, scp_fxx, gds_fxx, title=r'$F_{xx}$ absolute error', fname=pltdir+'errorfxx')

    time_ratio = scp_time[0] / gds_time[0]
    print(f"The series expansion method implemented in Fortran is {time_ratio:.0f} "
          f"times faster than the integral expressions evaluated with Scipy.")
