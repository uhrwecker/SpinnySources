import numpy as np
import matplotlib.pyplot as pl

from util.redshift import g
from util.solid_angles import alpha, beta

def get_redshift_area(ds):
    a = []
    b = []
    gg = []

    for quad in ds.get_all():
        for dp in quad.values():
            a.append(alpha(dp))
            b.append(beta(dp))
            gg.append(g(dp))

    a = np.array(a)
    b = np.array(b)
    gg = np.array(gg)

    arr_inds = a.argsort()

    a = a[arr_inds]
    b = b[arr_inds]
    gg = gg[arr_inds]

    b_mean = (b[a == np.amax(a)] + b[a == np.amin(a)]) / 2

    idx = np.where([b > b_mean[0]])[-1]
    id2 = np.where([b < b_mean[0]])[-1]

    xup = a[idx]
    yup = b[idx]
    gup = gg[idx]

    xdown = a[id2]
    ydown = b[id2]
    gdown = gg[id2]


    #matrix = []
    #last_gup = np.nan
    #for y in yup[np.argsort(yup)]:
    #    row = []
    #    for x in xup:
    #        if y == yup[x == xup]:
    #            row.append(gup[x == xup][0])
    #            last_gup = row[-1]
    #        else:
    #            row.append(np.nan)
    #    matrix.append(row)

    #from scipy.interpolate import interp2d
    #z = interp2d(xup, yup, gup)
    #z = z(np.linspace(np.amin(xup), np.amax(xup), 106), np.linspace(np.amin(yup), np.amax(yup), 106))
    #pl.pcolormesh(xup, yup[np.argsort(yup)], z)

    #pl.scatter(xup, yup, c=gup)
    #pl.scatter(xdown, ydown, c=gdown)

    pl.scatter(a, b, c=gg)
    pl.xlim(np.min(a)-0.1, np.max(a)+0.1)
    pl.ylim(np.min(b)-0.1, np.max(b)+0.1)
    pl.grid()
    pl.xlabel(r'$\alpha$')
    pl.ylabel(r'$\beta$')
    pl.show()

