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

    N1 = round(N * 0.30)
    N2 = N - N1

    N2a = round(N2 / 2)
    N2b = N2 - N2a

    # initial state
    system.gamma_list = [0.0] * system.num_slices

    # material properties
    slice_properties_1 = SliceProperties()
    slice_properties_1.ds_up = 80.0  # (eV^-1 nm^-3)
    slice_properties_1.ds_dn = 20.0  # (eV^-1 nm^-3)
    slice_properties_1.tau = 10000.0  # (fs)
    slice_properties_1.lifetime_up = 10000.0  # (fs)
    slice_properties_1.lifetime_dn = 10000.0  # (fs)

    slice_properties_2a = SliceProperties()
    slice_properties_2a.ds_up = 30.0  # (eV^-1 nm^-3)
    slice_properties_2a.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_2a.tau = 10000.0  # (fs)
    slice_properties_2a.lifetime_up = 10000.0  # (fs)
    slice_properties_2a.lifetime_dn = 10000.0  # (fs)

    slice_properties_2b = SliceProperties()
    slice_properties_2b.ds_up = 30.0  # (eV^-1 nm^-3)
    slice_properties_2b.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_2b.tau = 10000.0  # (fs)
    slice_properties_2b.lifetime_up = 10.0  # (fs)
    slice_properties_2b.lifetime_dn = 10.0  # (fs)

    plane_properties_1 = PlaneProperties()

    # case 1
    plane_properties_1.alpha_up = 8.0  # (eV^-1 nm^-1 fs^-1)
    plane_properties_1.alpha_dn = 2.0  # (eV^-1 nm^-1 fs^-1)

    # case 2
    #plane_properties_1.alpha_up = 9.0  # (eV^-1 nm^-1 fs^-1)
    #plane_properties_1.alpha_dn = 1.0  # (eV^-1 nm^-1 fs^-1)

    # case 3
    #plane_properties_1.alpha_up = 7.0  # (eV^-1 nm^-1 fs^-1)
    #plane_properties_1.alpha_dn = 3.0  # (eV^-1 nm^-1 fs^-1)

    plane_properties_2 = PlaneProperties()
    plane_properties_2.alpha_up = 0.4  # (eV^-1 nm^-1 fs^-1)
    plane_properties_2.alpha_dn = 0.4  # (eV^-1 nm^-1 fs^-1)

    system.slice_property_list = [slice_properties_1] * N1
    system.slice_property_list.extend([slice_properties_2a] * N2a)
    system.slice_property_list.extend([slice_properties_2b] * N2b)

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

    return system
