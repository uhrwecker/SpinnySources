import numpy as np


def geod(w, t):
    """
    Defines the differential equations for lightlike geodesics.

    Arguments:
        w :  vector of the state variables:
                  w = [t,t',r,r',theta,theta',phi,phi']
        t :  time
        m :  mass parameter
    """
    m = 1
    t, td, r, rd, th, thd, phi, phid = w

    alpha = 1 - 2 * m / r
    f = [td,
         - 2 * m * td * rd / (r ** 2 * alpha),
         rd,
         -m * td ** 2 / r ** 2 * alpha + m / r ** 2 * rd ** 2 / alpha + r * alpha * thd ** 2 + alpha * r * np.sin(
             th) ** 2 * phid**2,
         thd,
         -2 / r * rd * thd + np.cos(th) * np.sin(th) * phid ** 2,
         phid,
         - 2 / r * rd * phid - 2 / np.tan(th) * thd * phid ]

    return f