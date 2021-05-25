import numpy as np

from util.redshift import g
import matplotlib.pyplot as pl
from util.plot import UP_COLORS, DOWN_COLORS


def plot(ds, ax1, ax2, i_d=0):
    phi0 = []
    phi1 = []
    gup = []
    gdown = []

    for quad in ds.get_all_up():
        for dp in quad.values():
            phi0.append(dp.phi0)
            gup.append(g(dp))

    for quad in ds.get_all_down():
        for dp in quad.values():
            phi1.append(dp.phi0)
            gdown.append(g(dp))

    ax1.plot(phi0, gup, label=f's = {ds.s}, primary')

    ax1.set_xlim(0, np.pi * 2)
    ax1.set_ylabel('g')
    ax1.axhline(1, alpha=0.7, lw=1, color='black')
    ax1.legend()
    ax1.grid()

    ax2.plot(phi1, gdown, label=f's = {ds.s}, secondary')

    ax2.set_xlim(0, np.pi * 2)
    ax2.set_xlabel(r'$\varphi_0$')
    ax2.set_ylabel('g')
    ax2.axhline(1, alpha=0.7, lw=1, color='black')
    ax2.legend()
    ax2.grid()


def plot_one(ds):
    fig = pl.figure(figsize=(13, 8))
    ax1 = pl.gca()

    for d in ds:

        phi0 = []
        gup = []

        for quad in d.get_all_up():
            for dp in quad.values():
                phi0.append(dp.phi0)
                gup.append(g(dp))

        ax1.plot(phi0, gup, label=f's = {d.s}, primary')

    ticks = np.linspace(0, 2*np.pi, num=9, endpoint=True)
    labels = [r'{} $\pi$'.format(str(t / np.pi)[:4]) for t in ticks]

    pl.xticks(ticks, labels)

    ax1.set_xlim(0, np.pi * 2)
    ax1.set_ylabel('g')
    ax1.axhline(1, alpha=0.7, lw=1, color='black')
    ax1.legend()
    ax1.grid()

def main():
    from data.dataset import Dataset

    dss = []

    for s in [-1, 0, 1]:
        path = f'A:/Dokumente/Data/centre_geod/s{s}/'
        ds = Dataset(s, path)
        dss.append(ds)

    plot_one(dss)
    pl.show()

if __name__ == '__main__':
    main()