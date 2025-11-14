# This script creates points for testing gdeep.
import scipy as sc
import numpy as np
import numdifftools as nd


epsabs = 1e-10
epsrel = 1e-10


def fint(x, y):
    """Integral expression for F."""

    x2 = x**2
    ey = np.exp(-y)
    py = np.pi * ey
    dy = 2*ey

    phy0 = py*(sc.special.struve(0, x) + sc.special.y0(x))
    
    f = -dy * sc.integrate.quad(
        lambda t: np.exp(t) * (x2 + t**2)**-0.5,
        0.0,
        y,
        epsabs=epsabs,
        epsrel=epsrel,
    )[0] - phy0
    
    return f


def gfunc(p, q, k0):
    """Infinite-depth free-surface Green function."""
    
    z = p[2] + q[2]
    r = np.sqrt((p[0] - q[0])**2 + (p[1] - q[1])**2)
    rpq = np.sqrt(r**2 + (p[2] - q[2])**2)
    rpd = np.sqrt(r**2 + z**2)

    x = k0*r
    y = -k0*z
    
    g = 1/rpq + 1/rpd + k0*fint(x, y) - 1j*2*np.pi*k0*np.exp(-y)*sc.special.j0(x)

    return g


def gdeep(p, q, k0):
    """Green function and derivatives."""

    g_re = lambda f: gfunc(f, q, k0).real
    gradg_re = nd.Gradient(g_re)
    hessg_re = nd.Hessian(g_re)
    
    g_im = lambda f: gfunc(f, q, k0).imag
    gradg_im = nd.Gradient(g_im)
    hessg_im = nd.Hessian(g_im)
    
    g = g_re(p) + 1j*g_im(p)
    gradg = gradg_re(p) + 1j*gradg_im(p)
    hessg = hessg_re(p) + 1j*hessg_im(p)

    print(f'p = [{p[0]}_wp, {p[1]}_wp, {p[2]}_wp]')
    print(f'q = [{q[0]}_wp, {q[1]}_wp, {q[2]}_wp]')
    print(f'k0 = {k0}_wp')
    print(f'g_py = ({g.real:.16e}_wp, {g.imag:.16e}_wp)')
    print(f'gradg_py = [' + ', &\n'.join(f'({gradg[i].real:.16e}_wp, {gradg[i].imag:.16e}_wp)' for i in range(3)) + ']')
    print(f'hessg_py(1,:) = [' + ', &\n'.join(f'({hessg[0, i].real:.16e}_wp, {hessg[0, i].imag:.16e}_wp)' for i in range(3)) + ']')
    print(f'hessg_py(2,:) = [' + ', &\n'.join(f'({hessg[1, i].real:.16e}_wp, {hessg[1, i].imag:.16e}_wp)' for i in range(3)) + ']')
    print(f'hessg_py(3,:) = [' + ', &\n'.join(f'({hessg[2, i].real:.16e}_wp, {hessg[2, i].imag:.16e}_wp)' for i in range(3)) + ']')


if __name__ == "__main__":
    print('\n' + 'Region 1'.join(2*[40*'-']))
    p = np.array([1.1, -0.1, -0.5])
    q = np.array([1.2, 2.7, -6.9])
    k0 = 1.5
    gdeep(p, q, k0)

    print('\n' + 'Region 2'.join(2*[40*'-']))
    p = np.array([-0.8, 1.5, -1.9])
    q = np.array([-3.7, 2.2, -2.7])
    k0 = 2.3
    gdeep(p, q, k0)

    print('\n' + 'Region 4'.join(2*[40*'-']))
    p = np.array([0.9, 0.3, -1.6])
    q = np.array([-4.1, 3.9, -0.4])
    k0 = 1.3
    gdeep(p, q, k0)

    print('\n' + 'Region 4'.join(2*[40*'-']))
    p = np.array([2.6, -1.4, -3.8])
    q = np.array([-2.1, -1.8, -0.2])
    k0 = 2.9
    gdeep(p, q, k0)
