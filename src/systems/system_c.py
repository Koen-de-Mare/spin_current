from src.system import *


def make_system_c(h, dt, substeps, electrons_per_packet):
    system = System()

    # discretization parameters
    system.dt = dt
    system.substeps = substeps
    system.electrons_per_packet = electrons_per_packet  # (nm^-2)

    N = round(100.0 / h)
    system.num_slices = N
    system.slice_length = 100.0 / N

    N1 = round(N * 0.15)
    N2 = N - N1

    # initial state
    system.gamma_list = [0.0] * system.num_slices

    # material properties
    slice_properties_1 = SliceProperties()
    slice_properties_1.ds_up = 80.0  # (eV^-1 nm^-3)
    slice_properties_1.ds_dn = 20.0  # (eV^-1 nm^-3)
    slice_properties_1.tau = 100  # (NOT fs)

    slice_properties_2 = SliceProperties()
    slice_properties_2.ds_up = 30.0  # (eV^-1 nm^-3)
    slice_properties_2.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_2.tau = 100  # (NOT fs)

    plane_properties_1 = PlaneProperties()
    plane_properties_1.alpha_up = 6.0  # (eV^-1 nm^-1 fs^-1)
    plane_properties_1.alpha_dn = 4.0  # (eV^-1 nm^-1 fs^-1)

    plane_properties_2 = PlaneProperties()
    plane_properties_2.alpha_up = 0.4  # (eV^-1 nm^-1 fs^-1)
    plane_properties_2.alpha_dn = 0.4  # (eV^-1 nm^-1 fs^-1)

    system.slice_property_list = [slice_properties_1] * N1
    system.slice_property_list.extend([slice_properties_2] * N2)

    system.plane_property_list = [plane_properties_1] * N1
    system.plane_property_list.extend([plane_properties_2] * (N2 - 1))

    # laser properties
    system.t0 = 10.0  # (fs)
    system.pulse_duration = 20.0  # (fs), t_fwhm = 2 * sqrt(ln(2)) * pulse_duration = 1.6651 * pulse_duration
    system.peak_power = 5.0  # (eV nm^-2 fs^-1)
    system.penetration_depth: float = 5.0  # (nm)
    system.photon_energy = 1.0  # (eV)

    # hot electron properties
    system.vf = 1.0  # (nm fs^-1)
    system.lifetime_up = 10.0  # (fs)
    system.lifetime_dn = 5.0  # (fs)

    return system
