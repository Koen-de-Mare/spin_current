import matplotlib.pyplot as plt


from src.system import *


def plot_spin_current(system: System,  state_list, z: float, hot_resolved: bool):
    n = round(z / system.slice_length)

    time_list = []
    spin_current_tot = []
    spin_current_hot = []
    spin_current_cold = []

    for i in range(len(state_list)):
        (t, gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up,
         mu0_hot_dn, j_hot_up, j_hot_dn, j_up, j_dn, j_spin) = state_list[i]

        time_list.append(t)
        spin_current_tot.append(j_spin[n])
        spin_current_hot.append(j_hot_up[n] - j_hot_dn[n])
        spin_current_cold.append(j_up[n] - j_dn[n])

    if hot_resolved:
        plt.figure()

        plt.subplot(311)
        plt.plot(time_list, spin_current_tot)
        plt.ylabel('j_spin tot (nm^-2 fs^-1)')
        plt.grid(True)

        plt.subplot(312)
        plt.plot(time_list, spin_current_hot)
        plt.ylabel('j_spin hot (nm^-2 fs^-1)')
        plt.grid(True)

        plt.subplot(313)
        plt.plot(time_list, spin_current_cold)
        plt.ylabel('j_spin cold (nm^-2 fs^-1)')
        plt.xlabel('t (fs)')
        plt.grid(True)

        title = 'spin current at z = {} nm'.format(z)
        plt.suptitle(title)

        plt.show()
    else:
        plt.figure()

        #title = 'spin current at z = {} nm'.format(z)
        #plt.title(title)

        plt.plot(time_list, spin_current_tot)
        plt.ylabel('j_spin (nm^-2 fs^-1)')
        plt.xlabel('t (fs)')
        plt.xlim([0, 300.0])
        plt.ylim([0.0, 1.2])
        plt.grid(True)

        plt.show()