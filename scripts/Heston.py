import ctypes
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt


lib = ctypes.CDLL("lib_Heston.dll")


lib.Heston.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_int, ctypes.c_int
]
lib.Heston.restype = None

def run_simulation(S0, v0, mu, kappa, theta, sigma, rho, T, N, seed=None):
    if seed is None:
        seed = 1234

    precios = np.zeros(N + 1, dtype=np.double)
    vol = np.zeros(N + 1, dtype=np.double)

    lib.Heston(
        S0, v0, mu, kappa, theta, sigma, rho, T,
        vol.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        precios.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        N, seed
    )
    return precios, vol


def worker(args):
    return run_simulation(*args)


def simular_en_paralelo(n_simulaciones, S0, v0, mu, kappa, theta, sigma, rho, T, N):
    args_list = [(S0, v0, mu, kappa, theta, sigma, rho, T, N, 1234 + i) for i in range(n_simulaciones)]

    with ProcessPoolExecutor() as executor:
        resultados = list(executor.map(worker, args_list))

    precios_todos = np.array([r[0] for r in resultados])
    vols_todos = np.array([r[1] for r in resultados])
    return precios_todos, vols_todos


if __name__ == "__main__":
    S0 = 100
    v0 = 0.2
    mu = 0.05
    kappa = 2.0
    theta = 0.04
    sigma = 0.3
    rho = -0.7
    T = 1.0
    N = 252
    n_simulaciones = 1000

    print("Running Monte Carlo simulations...")
    precios, vols = simular_en_paralelo(n_simulaciones, S0, v0, mu, kappa, theta, sigma, rho, T, N)
    print("Simulaciones completadas.")

    # Graficar algunas trayectorias
    for i in range(100):
        plt.plot(precios[i], alpha=0.5)
    plt.xlabel("Time (days)")
    plt.ylabel("Asset Price")
    plt.grid(True)
    plt.show()