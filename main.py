from data.dataset import Dataset
from data.trajectory import TrajectoryResults
from results import plot_com as pc
from results import plot_redshift as pr
from results import plot_redshift_dist as prd
from util.redshift import g
from util.solid_angles import theta, phi, alpha, beta, get_phiminmax, get_thetaminmax
from util.edge_array import get_edge_points
from util.integrate_area import bounded_area
from util.collision import check_if_inside

import matplotlib.pyplot as pl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mp
import numpy as np

def plot_constants_of_motion(dss):
    fig1, (ax1, ax2, ax3) = pl.subplots(3, 1, figsize=(15, 15), sharex=True)


    for n, ds in enumerate(dss):
        pc.plot_com(ds, ax1, ax2, ax3, n)

    pl.show()

def plot_redshift_from_ds(dss):
    fig, (ax1, ax2) = pl.subplots(2, 1, figsize=(13, 7), sharex=True, sharey=True)

    for n, ds in enumerate(dss):
        pr.plot(ds, ax1, ax2, n)

    pl.show()

def plot_solid_angle(ds, both=True, ax=None):
    if not ax:
        fig = pl.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)

    colors = ['blue', 'red', 'green', 'orange']
    i = 0

    if both:
        for ds0 in ds:
            for quad, quad2 in zip(ds0.get_all_up(), ds0.get_all_down()):
                phi2 = []
                th2 = []
                phi0 = []
                cos = []

                for dp in quad.values():
                    phi2.append(alpha(dp))
                    th2.append(beta(dp))

                    #if np.isclose(dp.phi0, np.pi * 2) or np.isclose(dp.phi0, np.pi) or np.isclose(dp.phi0, 1.5 * np.pi):
                    #    ax.scatter(alpha(dp), beta(dp), marker='x', c=colors[i])

                for dp in quad2.values():
                    phi0.append(alpha(dp))
                    cos.append(beta(dp))

                ax.plot(phi2, th2, c=colors[i])
                ax.plot(phi0, cos, ls='--', c=colors[i])
            i += 1


    else:
        j = 0
        pmin, pmax = get_phiminmax(8, 0.5 * np.pi, np.pi / 2, 0.2)
        tmin, tmax = get_thetaminmax(8, 0.5 * np.pi, np.pi / 2, 0.2)
        print(pmin, pmax, tmin, tmax)
        pmax -= np.pi

        for ds0 in ds:
            all = []
            gs = []
            for quad in ds0.get_all():
                phi2 = []
                th2 = []

                for dp in quad.values():
                    #if check_if_inside(dp, tmin, tmax, pmin, pmax).size:
                    #    print(j)
                    #    j += 1
                    #else:
                    if dp.r0 >= 8.:
                        phi2.append(alpha(dp))
                        th2.append(beta(dp))
                        all.append([alpha(dp), beta(dp)])
                        gs.append(g(dp))

                ax.scatter(phi2, th2, c='black', s=1)
                print(len(phi2))
                #ax.scatter(0.03457183975193295, 7.159903237573876, c=colors[i], marker='x')

            i += 1

    p = np.linspace(0, np.pi * 2, num=1000)
    #ax.plot(np.cos(p), np.sin(p), c='black')

    #ax.set_xlim(-11, 11)
    #ax.set_ylim(-11, 11)

    #ax.set_xlabel(r'$\alpha$')
    #ax.set_ylabel(r'$\beta$')
    #ax.grid()
    #ax.legend()

    return np.array(all), np.array(gs), ax

def atan2(y, x, tol=1e-4):
    sig = []

    for xe, ye in zip(x, y):
        if -tol < xe < tol and -tol < ye < tol:
            sig.append(np.sign(ye) * np.pi / 2)
        elif xe > tol:
            sig.append(np.arctan(ye / xe))
        elif xe < tol:
            sig.append(np.arctan(ye / xe) + np.sign(ye) * np.pi)

    return np.array(sig)

def main():
    #path = 'A:/Dokumente/Data/centre_geod/'
    path = 'A:/Dokumente/Data/sphere_geod/phi_05_pi/'
    from scipy.interpolate import griddata

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = pl.subplots(3, 3, figsize=(10, 10), sharex=True,
                                                                           sharey=True)
    fp = [(-1, 's-1/', ax1), (-0.75, 's-075/', ax2), (-0.5, 's-05/', ax3),
          (-0.25, 's-025/', ax4), (0, 's0/', ax5), (0.25, 's025/', ax6),
          (0.5, 's05/', ax7), (0.75, 's075/', ax8), (1, 's1/', ax9)]

    cmap = pl.cm.cool_r
    norm = mp.colors.Normalize(1.0893496402614364, 1.1694971704654886)

    ga = []
    ass = []
    bss = []
    for s, f, ax in fp:
        d0 = TrajectoryResults('A:/Dokumente/Data/centre_geod/'+f+'1/up_1.57079633')
        g0 = g(d0)
        ds0 = Dataset(s, path+f, False)

        data, gs, ax = plot_solid_angle([ds0], False, ax)
        amin = np.amin(data[:, 0])
        amax = np.amax(data[:, 0])
        bmin = np.amin(data[:, 1])
        bmax = np.amax(data[:, 1])

        #gs /= g0

        ga.append(np.amin(gs))
        ga.append(np.amax(gs))

        ass.append(amin)
        ass.append(amax)
        bss.append(bmin)
        bss.append(bmax)

        gridx, gridy = np.mgrid[amin:amax:1000j, bmin:bmax:1000j]

        res = griddata(data, gs, (gridx, gridy), method='linear')

        ax.scatter(amin, bmin, s=0, label=f's = {s}')
        rim = ax.imshow(res.T, extent=(amin, amax, bmin, bmax), cmap=cmap, norm=norm)

        if s == -0.5 or s == 0.25 or s == 1:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            fig.colorbar(rim, cax=cax, orientation='vertical')
        ax.legend()#, gridx.T, gridx.T, shading='gouraud')

    ga = np.array(ga)
    ass = np.array(ass)
    bss = np.array(bss)

    ax.set_xlim(np.amin(ass), np.amax(ass))
    ax.set_ylim(np.amin(bss), np.amax(bss))

    ax1.set_ylabel(r'$\beta$')
    ax4.set_ylabel(r'$\beta$')
    ax7.set_ylabel(r'$\beta$')

    ax7.set_xlabel(r'$\alpha$')
    ax8.set_xlabel(r'$\alpha$')
    ax9.set_xlabel(r'$\alpha$')

    #fig.colorbar(rim, cax=axcb)

    print(np.amin(ass), np.amax(ass))
    print(np.amin(bss), np.amax(bss))
    print(np.amin(ga), np.amax(ga))

    pl.show()

def new():

    amin = None
    amax = None
    bmin = None
    bmax = None

    pmin, pmax = get_phiminmax(8, 1.5*np.pi, np.pi/2, 0.2)
    tmin, tmax = get_thetaminmax(8, 1.5*np.pi, np.pi/2, 0.2)

    for quad in ds0.get_all():
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

    from scipy.integrate import odeint
    from geodesics.equations import geod

    sigma = np.linspace(0, 25, num=100000)

    r = 8
    t = np.pi/2

    E = 0.7897629504703588
    L = -0.01824821023802678
    Q = 30.037438098323914

    q = Q / E**2
    l = L / E

    dt = 1 / (1 - 2/r)
    dr = -np.sqrt(1 - (q + l**2)/r**2 * (1 - 2/r))
    dtheta = -np.sqrt(q - l**2 / np.tan(t)) / r**2
    dphi = l / (r**2 * np.sin(t)**2)

    res = odeint(geod, [0, dt, r, dr, t, dtheta, 3*np.pi/2, dphi], sigma, atol=1e-5, rtol=1e-5)

    fig = pl.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    r = res[:, 2]
    t = res[:, 4]
    p = res[:, 6]

    ax.plot(r * np.cos(p) * np.sin(t),
            r * np.sin(p) * np.sin(t),
            r * np.cos(t))

    ax.scatter(0, -8, 0)
    ax.scatter(0, 15*np.sin(1.25), 15*np.cos(1.25))

    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_zlim(-10, 10)

    pl.show()


if __name__ == '__main__':
    main()