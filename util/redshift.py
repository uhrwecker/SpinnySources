import numpy as np


def g(datapoint):
    rem = datapoint.r[0]
    robs = datapoint.r[datapoint.r > 15][0]

    left = datapoint.gamma * np.sqrt(1 - 2 / robs)
    right = np.sqrt(1 - 2/rem) - datapoint.Lph / (datapoint.Eph * rem) * datapoint.vphi

    return (left * right)**(-1)
