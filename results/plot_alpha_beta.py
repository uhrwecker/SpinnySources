import numpy as np
import matplotlib.pyplot as pl

from util.solid_angles import alpha, beta

def plot(ds):
    fig = pl.figure(figsize=(13, 8))
    ax = pl.gca()

    a = []
    b = []

    for ups in ds.get_all_up():
        a = []
        b = []
        ps = []
        for point in ups.values():
            a.append(alpha(point))
            b.append(beta(point))
            ps.append(point.phi0)
        ps = np.array(ps)
        pmin = str(np.amin(ps) / np.pi)[:4]
        pmax = str(np.amax(ps) / np.pi)[:4]
        ax.plot(a, b, label=r'$\varphi$ from ' + pmin + r'$\pi$ to ' + pmax + r'$\pi$')

    a = np.array(a)
    b = np.array(b)

    ax.legend()


    p = np.linspace(0, np.pi * 2, num=1000)
    ax.plot(np.cos(p), np.sin(p), c='black', alpha=0.7)

    ax.grid()
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\beta$')


def main():
    from data.dataset import Dataset

    s = -1

    path = f'A:/Dokumente/Data/centre_geod/s{s}/'
    ds = Dataset(s, path)

    plot(ds)
    pl.show()

if __name__ == '__main__':
    main()