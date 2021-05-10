import numpy as np


def theta(dp):
    robs = dp.r[dp.r > 15][0]
    drobs = dp.dr[dp.r > 15][0]

    return np.arccos(drobs / dp.Eph)


def phi(dp):
    robs = dp.r[dp.r > 15][0]
    dphi = dp.dphi[dp.r > 15][0]
    t = dp.theta[dp.r > 15][0]

    return np.arcsin(dphi * robs * np.sin(t) / (dp.Eph * np.sin(theta(dp))) * np.sqrt(1 - 2 / robs))


def alpha(dp):
    robs = dp.r_obs
    thobs = dp.theta_obs
    root = np.sqrt(dp.Eph**2 / (1 - 2/robs) - (dp.Qph + dp.Lph**2) / robs**2)

    return dp.Lph / (np.sin(thobs) * root)


def beta(dp):
    factor = 1
    robs = dp.r_obs
    thobs = dp.theta_obs
    root = dp.Eph ** 2 / (1 - 2 / robs) - (dp.Qph + dp.Lph ** 2) / robs ** 2

    if dp.phi[0] < np.pi and dp.phi[0] > 0:
        factor = -1

    return factor * np.sqrt(np.abs((dp.Qph - dp.Lph**2 * np.tan(thobs)**(-2)) / root))


def get_phi(x, y, tol):
    sig = []

    for xe, ye in zip(x, y):
        if -tol < xe < tol and -tol < ye < tol:
            sig.append(np.sign(ye) * np.pi / 2)
        elif xe > tol:
            sig.append(np.arctan(ye / xe))
        elif xe < tol:
            sig.append(np.arctan(ye / xe) + np.sign(ye) * np.pi)

    return np.array(sig)


def get_phiminmax(r0, phi0, theta0, dr, tol=1e-5):
    p = np.linspace(0, 2*np.pi, num=1000, endpoint=False)
    t0 = np.pi / 2

    psi_x = r0 * np.cos(phi0) * np.sin(theta0) + dr * np.cos(p) * np.sin(t0)
    psi_y = r0 * np.sin(phi0) * np.sin(theta0) + dr * np.sin(p) * np.sin(t0)

    phi = get_phi(psi_x, psi_y, tol=tol)

    return np.amin(phi), np.amax(phi)


def get_thetaminmax(r0, phi0, theta0, dr):
    p = 0
    t0 = np.linspace(0, np.pi, num=1000)

    psi_x = r0 * np.cos(phi0) * np.sin(theta0) + dr * np.cos(p) * np.sin(t0)
    psi_y = r0 * np.sin(phi0) * np.sin(theta0) + dr * np.sin(p) * np.sin(t0)
    psi_z = r0 * np.cos(theta0) + dr * np.cos(t0)

    rho = np.sqrt(psi_x ** 2 + psi_y ** 2 + psi_z ** 2)
    theta = np.arccos(psi_z / rho)

    return np.amin(theta), np.amax(theta)
