import matplotlib.pyplot as plt


from src.system import *


def plot_spin_current(system: System, z: float, state_list):
    n = round(z / system.slice_length)

    time_list = []
    spin_current_list = []
    for i in range(len(state_list)):
        (t, gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up, \
         mu0_hot_dn, j_hot_up, j_hot_dn, j_up, j_dn, j_spin) = state_list[i]

        time_list.append(t)
        spin_current_list.append(j_spin[n])

    plt.plot(time_list, spin_current_list)
    plt.ylabel('j_spin (nm^-2 fs^-1')
    plt.xlabel('z (nm)')
    plt.grid(True)
    plt.show()