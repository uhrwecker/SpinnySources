import numpy as np

def a(r, t, l, q):
    root = l**2 / (np.sin(t)**2) * ((1 - 2/r)**(-1) - (q + l**2)/r**2)**(-1)
    return np.sqrt(root)

def b(r, t, l, q):
    root = (q - l**2 / np.tan(t)**2) / ((1 - 2/r)**(-1) - (q + l**2)/r**2)
    return -np.sqrt(root)

def lamda(r, t, al, be):
    return -al * np.sin(t) / np.sqrt(1 - 2/r) * np.sqrt(r**2 / (be**2 + al**2 + r**2))

def qu(r, t, al, be):
    l = -al * np.sin(t) / np.sqrt(1 - 2/r) * np.sqrt(r**2 / (be**2 + al**2 + r**2))
    return l**2 * (be**2 / (al**2 * np.sin(t)**2) + 1 / np.tan(t)**2)

def lamda2(r, t, al, be):
    return np.sqrt(qu(r, t, al, be) * al**2 * np.sin(t)**2 / (be**2 + al**2 * np.cos(t)**2))



def main():
    r = 15
    t = 1.25

    e = 0.7892856749702534
    l = -0.029047722002062425 / e
    q = 15.18665145735856 / e**2

    print(a(r, t, l, q), b(r, t, l, q))
    print(lamda(r, t, al=a(r, t, l, q), be=b(r, t, l, q)), qu(r, t, al=a(r, t, l, q), be=b(r, t, l, q)))
    print(l, q)

if __name__ == '__main__':
    main()