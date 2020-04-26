from src.system import *


def make_system_a(h, dt, electrons_per_packet):
    system = System()

    # discretization parameters
    system.dt = dt
    system.electrons_per_packet = electrons_per_packet  # (nm^-2)

    N = round(30.0 / h / 2)
    system.num_slices = 2 * N
    system.slice_length = 15.0 / N  # (nm)

    # initial state
    system.gamma_list = [0.0] * system.num_slices

    # material properties
    slice_properties_1 = SliceProperties()
    slice_properties_1.ds_up = 70.0  # (eV^-1 nm^-3)
    slice_properties_1.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_1.tau = 100  # (NOT fs)

    slice_properties_2 = SliceProperties()
    slice_properties_2.ds_up = 30.0  # (eV^-1 nm^-3)
    slice_properties_2.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_2.tau = 100  # (NOT fs)

    plane_properties_1 = PlaneProperties()
    plane_properties_1.alpha_up = 12.0  # (eV^-1 nm^-1 fs^-1)
    plane_properties_1.alpha_dn = 2.0  # (eV^-1 nm^-1 fs^-1)

    plane_properties_2 = PlaneProperties()
    plane_properties_2.alpha_up = 0.4  # (eV^-1 nm^-1 fs^-1)
    plane_properties_2.alpha_dn = 0.4  # (eV^-1 nm^-1 fs^-1)

    system.slice_property_list = [slice_properties_1] * N
    system.slice_property_list.extend([slice_properties_2] * N)

    system.plane_property_list = [plane_properties_1] * N
    system.plane_property_list.extend([plane_properties_2] * (N - 1))

    # laser properties
    system.t0 = 10.0  # (fs)
    system.pulse_duration = 5.0  # (fs)
    system.peak_power = 5.0  # (eV nm^-2 fs^-1)
    system.penetration_depth: float = 5.0  # (nm)
    system.photon_energy = 1.0  # (eV)

    # hot electron properties
    system.vf = 1.0  # (nm fs^-1)
    system.lifetime_up = 20.0  # (fs)
    system.lifetime_dn = 10.0  # (fs)

    return system
