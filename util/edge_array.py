import numpy as np
import matplotlib.pyplot as pl

from util.solid_angles import alpha, beta


def get_edge_points(ds):
    a = []
    b = []

    for quad in ds.get_all():
        for dp in quad.values():
            a.append(alpha(dp))
            b.append(beta(dp))

    a = np.array(a)
    b = np.array(b)

    # sort alpha and beta arrays
    arr_inds = a.argsort()

    a = a[arr_inds]
    b = b[arr_inds]

    # create lists for every quadrant
    a_ul = [np.amin(a)]
    b_ul = [b[a == np.amin(a)][0]]

    a_ll = [np.amin(a).tolist()]
    b_ll = [b[a == np.amin(a)][0]]

    a_ur = [np.amax(a).tolist()]
    b_ur = [b[a == np.amax(a)][0]]

    a_lr = [np.amax(a).tolist()]
    b_lr = [b[a == np.amax(a)][0]]

    for aa, bb in zip(a, b):
        # upper left:
        if np.amin(a) < aa <= a[b == np.amax(b)]:
            if b[a == np.amin(a)] < bb <= np.amax(b) and bb > b_ul[-1]:
                a_ul.append(aa)
                b_ul.append(bb)

            elif np.amin(b) < bb <= b[a == np.amin(a)] and bb < b_ll[-1]:
                a_ll.append(aa)
                b_ll.append(bb)

    for aa, bb in zip(a[::-1], b[::-1]):
        if a[b == np.amax(b)] < aa <= np.amax(a):
            if b[a == np.amax(a)] < bb <= np.amax(b) and bb > b_ur[-1]:# and bb > 0:
                a_ur.append(aa)
                b_ur.append(bb)

            elif np.amin(b) < bb <= b[a == np.amax(a)] and bb < b_lr[-1]:#
                a_lr.append(aa)
                b_lr.append(bb)

    alphas = np.array(a_ul + a_ur[::-1] + a_lr + a_ll[::-1])
    betas = np.array(b_ul + b_ur[::-1] + b_lr + b_ll[::-1])

    return alphas, betas
