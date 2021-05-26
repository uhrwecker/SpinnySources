import numpy as np
import matplotlib.pyplot as pl

from orbit_vel.circular import get_com, vphi

def plot(ss, r0=8):
    pl.figure(figsize=(13, 8))

    v = []

    for s in ss:
        E, _, L = get_com(s, r0)
        v.append(vphi(r0, s, L, E))

    pl.plot(ss, v, label=f'emitter at r0 = {r0}')
    pl.xlim(-1, 1)
    pl.xlabel('s')
    pl.ylabel(r'$v_\varphi$')
    pl.grid()
    pl.legend()


def main():
    s = np.linspace(-1, 1, num=1000)
    plot(s)
    pl.show()

if __name__ == '__main__':
    main()