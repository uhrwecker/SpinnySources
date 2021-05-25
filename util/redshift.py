import numpy as np


def g(datapoint):
    rem = datapoint.r[0]
    robs = datapoint.r[datapoint.r > 15][0]

    left = datapoint.gamma * np.sqrt(1 - 2 / robs)
    right = np.sqrt(1 - 2/rem) - datapoint.Lph / (datapoint.Eph * rem) * datapoint.vphi

    return (left * right)**(-1)


def g_from_l(l, rem, robs, vphi):
    gamma = 1/np.sqrt(1 - vphi**2)

    left = gamma * np.sqrt(1 - 2 / robs)
    right = np.sqrt(1 - 2 / rem) - l / rem * vphi

    return (left * right)**(-1)