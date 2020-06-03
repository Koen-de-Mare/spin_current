import matplotlib.pyplot as plt


import math


from src.injected import injected


def total_transfer_depthresolved(system, state_list, z_min: float, z_max: float, omega: float):
    z_list = []
    r_list = []
    t0_list = []

    prev_t0 = 50.0

    n_min = round(z_min / system.slice_length)
    n_max = round(z_max / system.slice_length)

    for n in range(n_min, n_max + 1):
        z = n * system.slice_length
        z_list.append(z)

        r, theta = injected(system, state_list, z, omega)
        t0 = theta / omega

        period = 2 * math.pi / omega
        t0 = t0 + round((prev_t0 - t0) / period) * period
        prev_t0 = t0

        r_list.append(r)
        t0_list.append(t0)


    plt.figure()

    plt.subplot(211)
    plt.plot(z_list, r_list)
    plt.ylabel('mode excitation (nm^-2)')
    plt.grid(True)

    plt.subplot(212)
    plt.plot(z_list, t0_list)
    plt.ylabel('t0 (fs)')
    plt.xlabel('z (nm)')
    plt.grid(True)

    title = 'depth dependence of excitation efficiency of a mode with omega = {} fs^-1'.format(omega)
    plt.suptitle(title)

    plt.show()