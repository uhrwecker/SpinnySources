import numpy as np
from util.plot import UP_COLORS, DOWN_COLORS

def plot_com(ds, ax1, ax2, ax3, i_d=0):
    phi0_up = []
    phi0_down = []
    Eph_up = []
    Eph_down = []
    Lph_up = []
    Lph_down = []
    Qph_up = []
    Qph_down = []

    for ups in ds.get_all_up():
        for point in ups.values():
            phi0_up.append(point.phi0)
            Eph_up.append(point.Eph)
            Lph_up.append(point.Lph)
            Qph_up.append(point.Qph)

    #for downs in ds.get_all_down():
    #    for point in downs.values():
    #        phi0_down.append(point.phi0)
    #        Eph_down.append(point.Eph)
    #        Lph_down.append(point.Lph)
    #        Qph_down.append(point.Qph)

    UP_COLORS = ['red', 'green', 'blue']
    ax1.plot(phi0_up, Eph_up, c=UP_COLORS[i_d], label=f's = {ds.s}, primary')
    #ax1.plot(phi0_down, Eph_down, c=DOWN_COLORS[i_d], label=f's = {ds.s}, secondary')
    ax1.set_xlim(0, 2*np.pi)
    ticks = np.linspace(0, 2 * np.pi, num=9, endpoint=True)
    labels = [r'{} $\pi$'.format(str(t / np.pi)[:4]) for t in ticks]

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)
    ax3.set_xlabel(r'$\varphi_0$')
    ax1.set_ylabel(r'$E_{photon}$')
    ax1.grid()
    ax1.legend()

    ax2.plot(phi0_up, Lph_up, c=UP_COLORS[i_d], label=f's = {ds.s}, primary')
    #ax2.plot(phi0_down, Lph_down, c=DOWN_COLORS[i_d], label=f's = {ds.s}, secondary')
    ax2.set_ylabel(r'$L_{photon}$')
    ax2.grid()
    ax2.legend()

    ax3.plot(phi0_up, Qph_up, c=UP_COLORS[i_d], label=f's = {ds.s}, primary')
    #ax3.plot(phi0_down, Qph_down, c=DOWN_COLORS[i_d], label=f's = {ds.s}, secondary')
    ax3.set_ylabel(r'$Q_{photon}$')
    ax3.grid()
    ax3.legend()


def plot_lambda_q(ds, ax1, ax2, i_d):
    phi0_up = []
    Eph_up = []
    Lph_up = []
    Qph_up = []

    for ups in ds.get_all_up():
        for point in ups.values():
            phi0_up.append(point.phi0)
            Eph_up.append(point.Eph)
            Lph_up.append(point.Lph)
            Qph_up.append(point.Qph)

    Eph_up = np.array(Eph_up)
    Lph_up = np.array(Lph_up)
    Qph_up = np.array(Qph_up)

    UP_COLORS = ['red', 'green', 'blue']
    ax1.plot(phi0_up, Lph_up / Eph_up, c=UP_COLORS[i_d], label=f's = {ds.s}, primary')
    ax1.set_xlim(0, 2 * np.pi)

    ticks = np.linspace(0, 2 * np.pi, num=9, endpoint=True)
    labels = [r'{} $\pi$'.format(str(t / np.pi)[:4]) for t in ticks]

    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)
    ax2.set_xlabel(r'$\varphi_0$')
    ax1.set_ylabel(r'$\lambda_{photon}$')
    ax1.grid()
    ax1.legend()

    ax2.plot(phi0_up, Qph_up / Eph_up**2, c=UP_COLORS[i_d], label=f's = {ds.s}, primary')
    ax2.set_ylabel(r'$q_{photon}$')
    ax2.grid()
    ax2.legend()


def main():
    import matplotlib.pyplot as pl
    from data.dataset import Dataset

    #fig1, (ax1, ax2, ax3) = pl.subplots(3, 1, figsize=(15, 15), sharex=True)
    fig1, (ax1, ax2) = pl.subplots(2, 1, figsize=(15, 15), sharex=True)

    dss = []

    for n, s in enumerate([-1, 0, 1]):
        path = f'A:/Dokumente/Data/centre_geod/s{s}/'
        ds = Dataset(s, path)
        dss.append(ds)

        #plot_com(ds, ax1, ax2, ax3, n)
        plot_lambda_q(ds, ax1, ax2, n)

    pl.show()

if __name__ == '__main__':
    main()