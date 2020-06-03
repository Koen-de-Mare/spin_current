import matplotlib.pyplot as plt


import math


from src.system import System


def plot_gamma(system: System, state_list, t_goal: float, zmin: float, zmax: float, gamma_min: float, gamma_max: float):
    i_best: int = 0
    cost: float = math.inf

    for i in range(len(state_list)):
        (t, gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up,
         mu0_hot_dn, j_hot_up, j_hot_dn, j_up, j_dn, j_spin) = state_list[i]

        cost_new = math.fabs(t_goal - t)

        if cost_new < cost:
            i_best = i
            cost = cost_new

    (t, gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up,
     mu0_hot_dn, j_hot_up, j_hot_dn, j_up, j_dn, j_spin) = state_list[i_best]

    (ticks_slices, ticks_planes) = system.make_ticks()

    plt.figure()

    plt.plot(ticks_slices, gamma)
    plt.xlim([zmin, zmax])
    plt.ylim([gamma_min, gamma_max])
    plt.ylabel('mu_s (eV)')
    plt.xlabel('z (nm)')
    plt.grid(True)

    title = 'gamma at = {} fs'.format(t)
    #plt.suptitle(title)

    plt.show()


