import numpy as np
import scipy.optimize as op

from util.collision import check_if_inside
from util.solid_angles import alpha, beta, get_phiminmax, get_thetaminmax
from scipy.integrate import odeint
from geodesics.equations import geod

def _get_min_max_alpha_beta(ds):
    amin = None
    amax = None
    bmin = None
    bmax = None

    for quad in ds.get_all():
        for dp in quad.values():
            a = alpha(dp)
            b = beta(dp)

            if not amin:
                amin = a
                amax = a
                bmax = b
                bmin = b

            if a < amin:
                amin = a
            if a > amax:
                amax = a
            if b < bmin:
                bmin = b
            if b > bmax:
                bmax = b

    return amin, amax, bmin, bmax

def _calc_l_q_from_alpha_beta(a, b, ig=(0, 0)):
    def func(w):
        l, q = w

        r = 15
        t = 1.25

        alp = 1 - 2/r
        gamma = np.sqrt(1 - 2/r) / np.sin(t)

        #first = gamma**2 * l**2 - a**2 * (1 - (q + l**2) / r**2 * (1 - 2/r))
        first = l**2 * (a**2 / r**2 * alp + gamma**2) + q * a**2 / r**2 - a**2
        #second = b **2 * (1 - (q + l**2) / r**2) - (1 - 2/r) * (q - l**2 / np.tan(t)**2)
        second = l**2 * (alp / np.tan(t)**2 - b**2 / r**2) - q * (b**2 / r**2 + alp) + b**2

        return first, second

    res = op.fsolve(func, x0=ig, xtol=1e-15)

    return res[0], res[1]


def calculate_area_geodesics(ds, phi0=1.5*np.pi):
    pmin, pmax = get_phiminmax(8, phi0, np.pi / 2, 0.2)
    tmin, tmax = get_thetaminmax(8, phi0, np.pi / 2, 0.2)

    amin, amax, bmin, bmax = _get_min_max_alpha_beta(ds)
    print(pmin, pmax, tmin, tmax)

    amin = -0.2
    amax = 0.2

    a = np.linspace(amin, amax, num=15)
    b = np.linspace(bmin, bmax, num=15)

    from data.trajectory import TrajectoryResults
    path = 'A:/Dokumente/Data/sphere_geod/phi_pi/'

    print(-4.9136823128767375/0.5734132847983956, 2.6858430669519446/0.5734132847983956**2)

    ls = []
    qs = []

    lm = []
    qm = []

    ll = []
    qq = []

    ig = (-10, 45)

    r = 15
    t = 1.25
    sigma = np.linspace(0, 35, num=10000)

    for aa in a:
        for bb in b:
            l, q = _calc_l_q_from_alpha_beta(aa, bb, ig=ig)

            dt = 1 / (1 - 2 / r)
            dr = -np.sqrt(1 - (q + l ** 2) / r ** 2 * (1 - 2 / r))
            dtheta = -np.sqrt(np.abs((q - l ** 2 / np.tan(t)) / r**4))
            dphi = l / (r ** 2 * np.sin(t) ** 2)

            res = odeint(geod, [0, dt, r, dr, t, dtheta, np.pi / 2, dphi], sigma, atol=1e-7, rtol=1e-7)

            #print(check_if_inside(res, tmin, tmax, pmin, pmax))

            if check_if_inside(res, tmin, tmax, pmin, pmax).size:
                ls.append(aa)
                qs.append(bb)
                #print(l, q)
            else:
                print(aa, bb)
                print(l, q)
                lm.append(aa)
                qm.append(bb)

                ll.append(l)
                qq.append(q)


        print(aa)

    #from scipy.integrate import odeint
    #from geodesics.equations import geod
    import matplotlib.pyplot as pl

    pl.scatter(ls, qs)
    pl.scatter(lm, qm)
    pl.show()

    sigma = np.linspace(0, 35, num=10000)

    fig = pl.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')
    q = 47.85229411526704
    l = 0.

    r = 15
    t = 1.25

    dt = 1 / (1 - 2 / r)
    dr = -np.sqrt(1 - (q + l ** 2) / r ** 2 * (1 - 2 / r))
    dtheta = -np.sqrt(np.abs((q - l ** 2 / np.tan(t)) / r**4))
    dphi = l / (r ** 2 * np.sin(t) ** 2)

    res = odeint(geod, [0, dt, r, dr, t, dtheta, np.pi/2, dphi], sigma, atol=1e-7, rtol=1e-7)

    r = res[:, 2]
    t = res[:, 4]
    p = res[:, 6]
    ax.plot(r * np.cos(p) * np.sin(t),
            r * np.sin(p) * np.sin(t),
            r * np.cos(t))

    ax.scatter(0, 15*np.sin(1.25), 15*np.cos(1.25))
    ax.scatter(0, -8, 0)

    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]
    x = 0.2 * np.cos(u) * np.sin(v)
    y = - 8 + 0.2 * np.sin(u) * np.sin(v)
    z = 0.2 * np.cos(v)
    ax.plot_wireframe(x, y, z, color="r")
    #psi_x = r0 * np.cos(phi0) * np.sin(theta0) + dr * np.cos(p) * np.sin(t0)
    #psi_y = r0 * np.sin(phi0) * np.sin(theta0) + dr * np.sin(p) * np.sin(t0)
    #psi_z = r0 * np.cos(theta0) + dr * np.cos(t0)

    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_zlim(-10, 10)

    pl.show()

if __name__ == '__main__':
    from data.dataset import Dataset
    path = 'A:/Dokumente/Data/sphere_geod/phi_3-2_pi/'

    ds0 = Dataset(0, path + 's0/', False)

    calculate_area_geodesics(ds0)