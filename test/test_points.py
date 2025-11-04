import os
import scipy as sc
import numpy as np
from time import perf_counter


epsabs = 1e-10
epsrel = 1e-10
spdir = 'scipy/'
os.makedirs(spdir, exist_ok=True)


def fint(x, y):
    """Integral expressions for F, Fx and Fxx."""

    x2 = x**2
    ey = np.exp(-y)
    py = np.pi * ey
    dy = 2*ey

    phy0 = py*(sc.special.struve(0, x) + sc.special.y0(x))
    phy1 = py*(sc.special.struve(1, x) + sc.special.y1(x))
    phy2 = py*(sc.special.struve(2, x) + sc.special.yn(2, x))
    
    f = -dy * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-0.5,
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] - phy0
    
    fx = dy*x * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-1.5,
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] - dy + phy1
    
    fxx = dy * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-2.5 * (t**2 - 2*x2),
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] + ey*x/3.0 + 0.5*(phy0 - phy2)
    
    return f, fx, fxx


def set_test_points(xlims=[0.005, 40.0], ylims=[0.005, 40.0], nx=500, ny=500):
    x1, x2 = xlims
    y1, y2 = ylims
    X, Y = np.mgrid[x1:x2:nx*1j, y1:y2:ny*1j]

    F = np.zeros_like(X)
    Fx = np.zeros_like(X)
    Fxx = np.zeros_like(X)

    for i, x in enumerate(X[:, 0]):
        for j, y in enumerate(Y[0, :]):
            f, fx, fxx = fint(x, y)
            F[i, j] = f
            Fx[i, j] = fx
            Fxx[i, j] = fxx

    np.save(spdir+'x.npy', X[:, 0])
    np.save(spdir+'y.npy', Y[0, :])
    np.save(spdir+'f.npy', F)
    np.save(spdir+'fx.npy', Fx)
    np.save(spdir+'fxx.npy', Fxx)


if __name__ == '__main__':
    start_time = perf_counter()
    
    set_test_points()
    
    elapsed_time = np.array([perf_counter() - start_time])
    np.save(spdir+'time.npy', elapsed_time)
    print(f'Elapsed time = {elapsed_time} seconds')
