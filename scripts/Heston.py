import ctypes
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
import os

dll_path = os.path.join(os.path.dirname(__file__), "lib_Heston.dll")
lib = ctypes.CDLL(dll_path)


lib.Heston.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_int, ctypes.c_int
]
lib.Heston.restype = None

dll_path = os.path.join(os.path.dirname(__file__), "HestonFFT_lib.dll")
heston_fft = ctypes.CDLL(dll_path)

heston_fft.CallHestonFFT.argtypes = [
    ctypes.c_double,  # Spot
    ctypes.c_double,  # Maturity
    ctypes.c_double,  # kappa
    ctypes.c_double,  # rho
    ctypes.c_double,  # volvol
    ctypes.c_double,  # theta
    ctypes.c_double,  # var0
    ctypes.c_double,  # rate
    ctypes.c_double,  # div
    ctypes.POINTER(ctypes.c_double),  # Strikes array
    ctypes.c_int,                     # numStrikes
    ctypes.POINTER(ctypes.c_double)  # CallPrices array (output)
]
heston_fft.CallHestonFFT.restype = None

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

def heston_call_prices_smile(
    S, T, r, v0, kappa, theta, sigma, rho,
    K_min=70, K_max=150, num=50, optType='Call', delta=0
):
    strikes = np.linspace(K_min, K_max, num)
    prices = []

    for K in strikes:
        price = lib.computeOptionPrice(S, K, r, delta, T, theta, kappa,
                                   sigma, rho, v0, optType.encode('utf-8'))
        prices.append(price)
    return prices, strikes


def heston_fft_prices(S, T, kappa, rho, sigma,
                      theta, v0, rate, div=0,
                      strikes=None, num_points=5000, strike_range=(80, 120)):
    # Generar strikes si no se proporcionan
    if strikes is None:
        strikes = np.linspace(strike_range[0], strike_range[1], num_points)
    else:
        strikes = np.array(strikes, dtype=np.float64)
    
    num_strikes = len(strikes)
    call_prices = np.zeros(num_strikes, dtype=np.float64)

    strikes_ptr = strikes.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    prices_ptr = call_prices.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    heston_fft.CallHestonFFT(
        S, T, kappa, rho, sigma,
        theta, v0, rate, div,
        strikes_ptr, num_strikes, prices_ptr
    )

    return strikes, call_prices


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