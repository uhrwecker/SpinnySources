'''Routine to calculate the area that is enclosed by an array'''

import numpy as np
import scipy.integrate as inter

def bounded_area(x, y):
    y += np.abs(np.amin(y))

    up_y = []
    up_x = []
    low_y = []
    low_x = []

    upper = True

    for xx, yy in zip(x, y):
        if upper:
            up_y.append(yy)
            up_x.append(xx)
        else:
            low_y.append(yy)
            low_x.append(xx)

        if xx == np.amax(x):
            upper = False

    up = inter.simpson(up_y, up_x)
    low = inter.simpson(low_y, low_x)

    return up - np.abs(low)

