import numpy as np
fr = -1
ft = -1

def geo(w, t, e=0, l=0, q=0):
    t, r, theta, phi = w

    alpha = 1 - 2 / r

    global fr, ft

    rroot = e**2 - (q + l ** 2) / r ** 2 * alpha

    if rroot < 8e-5:
        print(rroot)
        fr *= -1

    #if np.abs(e**2 - (q + l ** 2) / r ** 2 * alpha) < 1e-3:
    #    print(fr * np.sqrt(1 - (q + l ** 2) / r ** 2 * alpha))
    #    fr *= -1

    troot = q - l**2 * np.tan(theta)**(-2)
    print(troot)
    if troot < 5e-2:
        #troot = np.abs(troot)
        print(troot)
        ft *= -1

    #if 0 < ft * np.sqrt(q - l**2 * np.tan(theta)**(-2)) / r**2 < 1e-04:
    #    print(np.sqrt(q - l**2 * np.tan(theta)**(-2)) / r**2)
    #    ft *= -1

    f = [
        e / alpha,
        fr * np.sqrt(e**2 - (q + l**2) / r**2 * alpha),
        ft * np.sqrt(troot) / r**2,
        l / (r**2 * np.sin(theta)**2)
    ]

    return f


def geo2(w, t):
    q = 8.168553024175978
    global l

    r, theta, phi = w

    global fr, ft

    if -2e-2 < 2 * l * (1 - 2/r) / r**2 * np.sqrt(1 - (q + l**2) / r**2 * (1 - 2/r)) < 0:
        fr *= -1

    f = [
        fr * 2 * l * (1 - 2/r) / r**2 * np.sqrt(1 - (q + l**2) / r**2 * (1 - 2/r)),
        ft * 2 * l / (r**2 * np.tan(theta)**2 * np.sqrt(q - l**2 / np.tan(theta)**2)),
        1 / (r**2 * np.sin(theta)**2)
    ]
    print(f)
    return f


def main():
    from scipy.integrate import odeint
    import matplotlib.pyplot as pl
    from mpl_toolkits.mplot3d import Axes3D

    fig = pl.figure()
    ax = fig.add_subplot(111, projection='3d')

    for f in [1]:
        sigma = np.linspace(0, 20, num=100000) * f
        res2 = odeint(geo, [0, 8, np.pi/2, np.pi], sigma, atol=1e-07, rtol=1e-07)
        #print(res2[:, 2])

        r = res2[:, 1]
        t = res2[:, 2]
        p = res2[:, 3]

        ax.plot(r * np.cos(p) * np.sin(t),
                r * np.sin(p) * np.sin(t),
                r * np.cos(t))

    sigma = np.linspace(0, 20, num=100000)
    res = odeint(geo2, [8, np.pi/2, np.pi], sigma, atol=1e-07, rtol=1e-07)

    r = res2[:, 1] + res[:, 0]-8
    t = res2[:, 2] + res[:, 1]-np.pi/2
    p = res2[:, 3] + res[:, 2]-np.pi

    ax.plot(r * np.cos(p) * np.sin(t),
            r * np.sin(p) * np.sin(t),
            r * np.cos(t))

    ax.scatter(15 * np.cos(np.pi / 2) * np.sin(1.25),
               15 * np.sin(np.pi / 2) * np.sin(1.25),
               15 * np.cos(1.25))

    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_zlim(-15, 15)

    pl.show()