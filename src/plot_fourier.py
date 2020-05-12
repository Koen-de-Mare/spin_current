import matplotlib.pyplot as plt


import math


from src.system import System
from src.injected import injected


def plot_fourier(system: System, state_list, z: float, omega_max: float):
    #  omega_max in (fs^-1)

    omega_list = []  # (fs^-1)
    r_list = []      # (nm^-2)
    t0_list = []     # (fs)

    prev_t0 = 50.0  # (fs)

    num_steps = 1000
    for i in range(1, num_steps):
        omega = omega_max * i / num_steps
        omega_list.append(omega)

        r, theta = injected(system, state_list, z, omega)
        t0 = theta / omega  # (fs)

        period = 2 * math.pi / omega
        t0 = t0 + round((prev_t0 - t0) / period) * period
        prev_t0 = t0

        r_list.append(r)
        t0_list.append(t0)

    plt.figure()

    plt.subplot(211)
    plt.plot(omega_list, r_list)
    plt.ylabel('r (nm^-2)')
    plt.grid(True)

    plt.subplot(212)
    plt.plot(omega_list, t0_list)
    plt.ylabel('t0 (fs)')
    plt.xlabel('omega (fs^-1)')
    plt.grid(True)

    title = 'fourier transformed spin current at z = {} nm'.format(z)
    plt.suptitle(title)

    plt.show()
