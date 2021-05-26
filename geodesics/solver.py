import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as pl

from geodesics.equations import geod


def solve(t, r, theta, phi, l, q, fr=-1, ft=1):
    sigma = np.linspace(0, 25, num=10000)

    dt = (1 - 2/r)**(-1)
    dr = fr * np.sqrt(1 - (q + l**2)/r**2 * (1 - 2/r))
    dtheta = ft * np.sqrt(np.abs(q - l**2 / np.tan(theta)**2)) / r**2
    dphi = l / (r**2 * np.sin(theta)**2)

    res = odeint(geod, [t, dt, r, dr, theta, dtheta, phi, dphi], sigma, atol=1e-7, rtol=1e-7)

    return res

    from mpl_toolkits.mplot3d import Axes3D

    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')

    r = res[:, 2]
    t = res[:, 4]
    p = res[:, 6]
    ax.plot(r * np.cos(p) * np.sin(t),
            r * np.sin(p) * np.sin(t),
            r * np.cos(t))

    ax.scatter(0, 15 * np.sin(1.25), 15 * np.cos(1.25))


    u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:10j]

    r0 = 8
    p0 = 0
    t0 = np.pi / 2

    x = 8 * np.cos(p0) * np.sin(t0) + 0.2 * np.cos(u) * np.sin(v)
    y = 8 * np.sin(p0) * np.sin(t0) + 0.2 * np.sin(u) * np.sin(v)
    z = 8 * np.cos(t0) + 0.2 * np.cos(v)

    ax.plot_wireframe(x, y, z, color="r")

    ax.scatter(8 * np.cos(p0), 8 * np.sin(p0), 0)

    #ax.set_xlim(-1, 1)
    #ax.set_ylim(-9, -7)
    #ax.set_zlim(-1, 1)

    pl.show()


#e = 0.7892856749702534
#l = -0.029047722002062425 / e
#q = 15.18665145735856 / e**2

#l = -0.1916548077951617
#q = 30.152133903045183

#l = -0.23977054367329073
#q = 29.744250501787644

#l = 0.007644096289870729
#q = 46.39294016841925

#l = -8.540143666368598
#q = 8.052409102365688
l = 8.57137871443015
q = 8.160479792842649

solve(0, 15, 1.25, np.pi/2, -l, q, ft=-1)
#1.5458015331759765 1.5957911204138167