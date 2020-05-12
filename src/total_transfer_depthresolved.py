import matplotlib.pyplot as plt

from src.injected import injected

def total_transfer_depthresolved(system, state_list, z_min: float, z_max: float, omega: float):
    z_list = []
    j_spin_list = []

    n_min = round(z_min / system.slice_length)
    n_max = round(z_max / system.slice_length)

    for n in range(n_min, n_max + 1):
        z = n * system.slice_length
        z_list.append(z)

        r, theta = injected(system, state_list, z, omega)
        j_spin_list.append(r)
        
    plt.figure()
    plt.title('depth dependence of excitation efficiency of a mode with omega = {}'.format(omega))
    plt.plot(z_list, j_spin_list)
    plt.ylabel('mode excitation (nm^-2)')
    plt.xlabel('z (nm)')
    plt.show()