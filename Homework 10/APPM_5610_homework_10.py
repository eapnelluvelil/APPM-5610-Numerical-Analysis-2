from re import M
from scipy.fft import dstn, idstn
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def main():
    # f = lambda x, y: np.sin(np.pi*x)*np.sin(np.pi*y)
    # f = lambda x, y: np.sin(4*np.pi*x)*np.sin(np.pi*y)
    f = lambda x, y: -np.sin(np.pi*x)*np.sin(5*np.pi*y)

    # Discretize grid
    a = 0
    b = 1
    p = 7
    # To efficiently use the Cooley-Tukey algorithm, we use 2**p points in
    # either direction
    q = 2**p

    x = np.linspace(a, b, q)
    y = x
    X, Y = np.meshgrid(x, y)

    # Compute sine wavenumbers
    m = np.arange(1, q + 1)
    n = m
    M, N = np.meshgrid(m, n)

    # Compute u at gridpoints via the DST and IDST
    f_hat = dstn(f(X, Y), norm="ortho")
    u_hat = np.divide(-f_hat, (np.pi)**2*(M**2 + N**2))
    u = idstn(u_hat, norm="ortho")

    # Plot solution
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(X, Y, u, cmap="coolwarm",
                           linewidth=0, antialiased=True)
    fig.colorbar(surf, shrink=0.4, aspect=5)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_zlabel("$u ( x, y )$")
    plt.show()

main()