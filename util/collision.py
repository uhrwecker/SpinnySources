import numpy as np


def check_if_inside(dp, tmin, tmax, pmin, pmax):
    rmin = 7.8
    rmax = 8.2

    if type(dp) == np.ndarray:
        dpr = dp[:, 2]
        dptheta = dp[:, 4]
        dpphi = dp[:, 6]
    else:
        dpr = dp.r
        dptheta = dp.theta
        dpphi = dp.phi

    r = dpr[dpr < rmax]
    r2 = r[r > rmin]

    idx1 = np.array([np.where(dpr == r2[n])[0][0] for n in range(len(r2))])

    t = dptheta[dptheta > tmin]
    t2 = t[t < tmax]

    idx2 = np.array([np.where(dptheta == t2[n])[0][0] for n in range(len(t2))])

    p = dpphi[dpphi > pmin]
    p2 = p[p < pmax]

    idx3 = np.array([np.where(dpphi == p2[n])[0][0] for n in range(len(p2))])

    a = np.intersect1d(idx1, idx2)
    b = np.intersect1d(idx2, idx3)
    c = np.intersect1d(idx3, idx1)

    if a.size and b.size and c.size:
        print(r[np.intersect1d(a, idx3)])
        print(t[np.intersect1d(a, idx3)])
        print(p[np.intersect1d(a, idx3)])
        return np.intersect1d(a, idx3)
    else:
        return np.array([])

#from data.trajectory import TrajectoryResults
#from util.solid_angles import get_phiminmax, get_thetaminmax##

#dp = TrajectoryResults('A:/Dokumente/Data/sphere_geod/phi_05_pi/s0/8.18578505_1.57623312_1.56343536')
#pmin, pmax = get_phiminmax(8, 0.5 * np.pi, np.pi / 2, 0.2)
#tmin, tmax = get_thetaminmax(8, 0.5 * np.pi, np.pi / 2, 0.2)

#print(tmin, tmax)
#print(pmin, (pmax+np.pi)%(2*np.pi))
#print(check_if_inside(dp, tmin, tmax, tmin, tmax))
