import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mp
from mpl_toolkits.axes_grid1 import make_axes_locatable

from results.find_geod_by_backwards import generate_hit_matrix
from data.trajectory import TrajectoryResults


def _set_ax_lim(ax, xmin, xmax, ymin, ymax):
    if np.abs(xmax - xmin) > np.abs(ymax - ymin):
        mean = (ymin + ymax) / 2
        dist = np.abs(xmax - xmin) / 2
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(mean - dist, mean + dist)
    else:
        mean = (xmin + xmax) / 2
        dist = np.abs(ymax - ymin) / 2
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(mean - dist, mean + dist)

def plot_matrix(ax, fig, mat, amin, amax, bmin, bmax, cmap, norm, s):
    rim = ax.imshow(mat, extent=(amin, amax, bmin, bmax), cmap=cmap, norm=norm)
    ax.scatter(amin, bmin, s=0, label=f's = {s}')

    if s == -0.5 or s == 0.25 or s == 1:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(rim, cax=cax, orientation='vertical')
    ax.legend()

    _set_ax_lim(ax, amin, amax, bmin, bmax)
    #ax.set_xlim(amin, amax)
    #ax.set_ylim(bmin, bmax+1)


def main():
    r = 15
    t = 1.25

    bmin = -0.3
    bmax = 0.3
    amax = 10.5
    amin = 9.8

    NUM = 70

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = pl.subplots(3, 3, figsize=(12, 12), sharex=True,
                                                                           sharey=True)

    #fig, ax = pl.subplots(1, 1, figsize=(10, 10))

    cmap = pl.cm.cool_r
    norm = mp.colors.Normalize(0.7197263059861253, 0.7813565661483185)

    fff = [(-1, 's-1/', ax1), (-0.75, 's-075/', ax2), (-0.5, 's-05/', ax3),
          (-0.25, 's-025/', ax4), (0, 's0/', ax5), (0.25, 's025/', ax6),
          (0.5, 's05/', ax7), (0.75, 's075/', ax8), (1, 's1/', ax9)]

    gs = []
    fs = []

    for s, fp, ax in fff:
        print(s)
        path = f'A:/Dokumente/Data/centre_geod/{fp}2/up_3.14159265'
        dp = TrajectoryResults(path)
        print(dp.Lph / dp.Eph)
        print(dp.Qph / dp.Eph**2)

        mat = generate_hit_matrix(dp.vphi, phi=dp.phi0, ft=-1, bmin=bmin, bmax=bmax, amin=amin, amax=amax, num=NUM)

        gs.append(np.amax(np.array(mat)[~np.isnan(mat)].flatten()))
        gs.append(np.amin(np.array(mat)[~np.isnan(mat)].flatten()))

        plot_matrix(ax, fig, mat, amin, amax, bmin, bmax, cmap, norm, s)

        flux = []
        for row in mat:
            idx = ~np.isnan(row)
            if True in idx:
                flux.append(np.trapz(np.array(row)[idx]**3, np.linspace(amin, amax, num=NUM)[idx]))
            else:
                flux.append(np.nan)

        idx = ~np.isnan(np.array(flux))
        fs.append(np.abs(np.trapz(np.array(flux)[idx], np.linspace(bmin, bmax, num=NUM)[idx])))

    gs = np.array(gs)

    print(np.amax(gs), np.amin(gs))

    pl.figure(figsize=(13, 6))
    pl.scatter([-1, -0.75, -0.5, -0.25, -0, 0.25, 0.5, 0.75, 1], fs)
    pl.grid()
    pl.xlabel('s')
    pl.ylabel('$F_{obs}$')

    pl.show()




if __name__ == '__main__':
    main()