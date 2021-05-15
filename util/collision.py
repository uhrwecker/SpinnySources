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
        #print(p[np.intersect1d(a, idx3)])
        #print(r[np.intersect1d(a, idx3)])
        #print(t[np.intersect1d(a, idx3)])
        return np.intersect1d(a, idx3)
    else:
        #print(dpr)#

        #print(a, b, c)
        return np.array([])
