import math


from src.system import System


def injected(system: System, state_list, z: float, omega: float) -> (float, float):
    # omega in (fs^-1)

    n = round(z / system.slice_length)

    a = 0.0  # (nm^-2)
    b = 0.0  # (nm^-2)

    for i in range(len(state_list)):
        (t, gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up,
         mu0_hot_dn, j_hot_up, j_hot_dn, j_up, j_dn, j_spin) = state_list[i]

        a += j_spin[n] * math.cos(t * omega) * system.dt
        b += j_spin[n] * math.sin(t * omega) * system.dt

    r = math.sqrt(a*a + b*b)  # (nm^-2)
    theta = math.atan2(b, a)  # (1)

    return r, theta
