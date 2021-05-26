import numpy as np
import time

from util.solid_angles import lamda, qu
from util.redshift import g_from_l
from geodesics.solver import solve


def check_if_hit(res, phi=np.pi/2):
    r = res[:, 2]
    t = res[:, 4]
    p = res[:, 6]

    x0 = 8 * np.cos(phi) * np.sin(np.pi/2)
    y0 = 8 * np.sin(phi) * np.sin(np.pi/2)
    z0 = 8 * np.cos(np.pi/2)

    x = r * np.cos(p) * np.sin(t)
    y = r * np.sin(p) * np.sin(t)
    z = r * np.cos(t)

    dist = np.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
    #print(x[np.logical_and(r < 8.2, r > 7.8)], x0)
    #print(y[np.logical_and(r < 8.2, r > 7.8)], y0)
    #print(z[np.logical_and(r < 8.2, r > 7.8)], z0)

    if dist[dist <= 0.2].size:
        return True
    else:
        return False


def generate_hit_matrix(vphi, phi=np.pi/2, fr=-1, ft=1,
                        bmin=-5.5, bmax=-4, amin=-0.6, amax=0.6, num=100):
    r0 = 15
    t0 = 1.25
    phi0 = np.pi/2

    mat = []

    ass = np.linspace(amin, amax, num=num)
    bss = np.linspace(bmax, bmin, num=num)

    for beta in bss:
        row = []
        for alpha in ass:
            l = lamda(r0, t0, alpha, beta)
            q = qu(r0, t0, alpha, beta)

            #if phi <= (np.pi+1e-2):
            #    l *= -1

            res = solve(0, r0, t0, phi0, -l, q, fr=fr, ft=ft)
            num = check_if_hit(res, phi)
            if num:
                row.append(g_from_l(l, 8, 15, vphi))
            else:
                row.append(np.nan)
        mat.append(row)

    return mat


def main():
    import matplotlib.pyplot as pl

    vphi = 0.40824829045486405

    mat = generate_hit_matrix(vphi, num=100)

    pl.imshow(mat, extent=(-0.6, 0.6, -5.5, -4.))
    pl.show()

if __name__ == '__main__':
    start = time.time()
    main()
    print(f'Took {time.time() - start}s.')