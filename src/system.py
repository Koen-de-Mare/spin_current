import matplotlib.pyplot as plt

import math
import random

electrons_per_packet = 0.001  # (nm^-2)
vf = 1.0  # (nm fs^-1)
penetration_depth = 5.0  # (nm)
photon_energy = 1.0  # (eV)

lifetime_up = 1000.0  # (fs)
lifetime_dn = 10.0  # (fs)


class HotElectronPacket:
    def __init__(self):
        self.is_up = True
        self.z = 0.0
        self.vz = 0.0


class SliceProperties:
    def __init__(self):
        self.ds_up = 0
        self.ds_dn = 0
        self.tau = 0


class PlaneProperties:
    def __init__(self):
        self.alpha_up = 0
        self.alpha_dn = 0


class System:
    def __init__(self):
        # state
        self.gamma_list: [float] = []
        self.t: float = 0.0
        self.hot_list = []

        # state dynamics
        self.j_hot: [float] = []
        self.j_up: [float] = []
        self.j_dn: [float] = []

        # properties
        self.num_slices = 0
        self.slice_length: float = 0.0  # (nm)
        self.slice_property_list: [SliceProperties] = []
        self.plane_property_list: [PlaneProperties] = []
        self.dt = 0.0  # (fs)

    def step(self):
        self.t += self.dt

        # transport of hot electrons -----------------------------------------------------------------------------------
        self.j_hot = [0.0] * (self.num_slices - 1)

        for i in range(len(self.hot_list)):
            hot_electron = self.hot_list[i]
            z0 = hot_electron.z
            hot_electron.z += hot_electron.vz * self.dt
            if hot_electron.z < 0.0:
                hot_electron.z = -hot_electron.z
                hot_electron.vz = -hot_electron.vz
            if hot_electron.z > self.slice_length * self.num_slices:
                hot_electron.z = 2 * self.slice_length * self.num_slices - hot_electron.z
                hot_electron.vz = -hot_electron.vz
            z1 = hot_electron.z

            zmin = min(z0, z1)
            jmin = math.ceil(zmin / self.slice_length - 1.0)
            zmax = max(z0, z1)
            jmax = math.ceil(zmax / self.slice_length - 1.0)

            assert(jmin >= 0)
            assert(jmax < self.num_slices)

            sign = 0
            if z1 > z0:
                sign = 1.0
            else:
                sign = -1.0

            for j in range(jmin, jmax):
                self.j_hot[j] += sign * electrons_per_packet

        # motion of thermal electrons ----------------------------------------------------------------------------------
        self.j_up = [0.0] * (self.num_slices - 1)
        self.j_dn = [0.0] * (self.num_slices - 1)

        for i in range(self.num_slices - 1):
            j_up_0_i = -self.plane_property_list[i].alpha_up * (
                    self.gamma_list[i+1] * self.slice_property_list[i+1].ds_dn /
                        (self.slice_property_list[i+1].ds_up + self.slice_property_list[i+1].ds_dn) -
                    self.gamma_list[ i ] * self.slice_property_list[ i ].ds_dn /
                        (self.slice_property_list[ i ].ds_up + self.slice_property_list[ i ].ds_dn)
                ) / self.slice_length
            j_dn_0_i = +self.plane_property_list[i].alpha_dn * (
                    self.gamma_list[i+1] * self.slice_property_list[i+1].ds_up /
                        (self.slice_property_list[i+1].ds_up + self.slice_property_list[i+1].ds_dn) -
                    self.gamma_list[ i ] * self.slice_property_list[ i ].ds_up /
                        (self.slice_property_list[ i ].ds_up + self.slice_property_list[ i ].ds_dn)
                ) / self.slice_length

            ee_i = \
                (j_up_0_i + j_dn_0_i + self.j_hot[i]) / \
                (self.plane_property_list[i].alpha_up + self.plane_property_list[i].alpha_dn)

            self.j_up[i] = j_up_0_i - self.plane_property_list[i].alpha_up * ee_i
            self.j_dn[i] = j_dn_0_i - self.plane_property_list[i].alpha_dn * ee_i

        # time derivative of gamma
        gamma_dot = [0.0] * self.num_slices

        for i in range(self.num_slices):
            gamma_dot[i] -= self.gamma_list[i] * (
                1.0 / self.slice_property_list[i].ds_up + 1.0 / self.slice_property_list[i].ds_dn
            ) / self.slice_property_list[i].tau

        for i in range(self.num_slices - 1):
            gamma_dot[i] += \
                (
                    -self.j_up[i] / self.slice_property_list[i].ds_up +
                     self.j_dn[i] / self.slice_property_list[i].ds_dn
                ) / self.slice_length
            gamma_dot[i + 1] += \
                (
                    self.j_up[i] / self.slice_property_list[i+1].ds_up -
                    self.j_dn[i] / self.slice_property_list[i+1].ds_dn
                ) / self.slice_length

        for i in range(0, self.num_slices):
            self.gamma_list[i] += self.dt * gamma_dot[i]

        # excitation and decay of hot electrons ------------------------------------------------------------------------

        # accumulators for net excitation of thermal electrons
        excited_packets_up = [0] * self.num_slices
        excited_packets_dn = [0] * self.num_slices

        new_hot_list = []

        # decay
        for i in range(len(self.hot_list)):
            lifetime = 0
            if self.hot_list[i].is_up:
                lifetime = lifetime_up
            else:
                lifetime = lifetime_dn

            if random.random() < math.exp(-self.dt / lifetime):
                # keep packet
                new_hot_list.append(self.hot_list[i])
            else:
                # remove packet
                slice_index = math.floor(self.hot_list[i].z / self.slice_length)
                if self.hot_list[i].is_up:
                    excited_packets_up[slice_index] -= 1
                else:
                    excited_packets_dn[slice_index] -= 1

        # excitation
        fluence = 10.0 * math.exp(-self.t / 10.0)  # (eV nm^-2 fs^-1)

        for i in range(self.num_slices):
            ds_tot_i = self.slice_property_list[i].ds_up + self.slice_property_list[i].ds_dn

            # (eV nm^-2 fs^-1)
            power_i = \
                -fluence * \
                math.exp(-self.slice_length * i / penetration_depth) * \
                (math.exp(-self.slice_length / penetration_depth) - 1.0)

            excited_packets_i = round(power_i * self.dt / (photon_energy * electrons_per_packet))
            excited_packets_i_up = round(
                excited_packets_i * self.slice_property_list[i].ds_up / ds_tot_i
            )
            excited_packets_i_dn = round(
                excited_packets_i * self.slice_property_list[i].ds_dn / ds_tot_i
            )

            excited_packets_up[i] += excited_packets_i_up
            excited_packets_dn[i] += excited_packets_i_dn

            for j in range(excited_packets_i_up + excited_packets_i_dn):
                new_packet = HotElectronPacket()

                new_packet.is_up = j < excited_packets_i_up
                new_packet.z = self.slice_length * (i + random.random())
                new_packet.vz = vf * (2.0 * random.random() - 1.0)

                new_hot_list.append(new_packet)

        self.hot_list = new_hot_list

        # apply the change in number of thermal electrons caused by excitation and decay
        for i in range(self.num_slices):
            self.gamma_list[i] += \
                electrons_per_packet * (
                    -excited_packets_up[i] / self.slice_property_list[i].ds_up +
                    excited_packets_dn[i] / self.slice_property_list[i].ds_dn
                )

    def make_data(self):
        # hot electron density
        hot_up = [0] * self.num_slices
        hot_dn = [0] * self.num_slices

        for packet in self.hot_list:
            slice_index = math.floor(packet.z / self.slice_length)
            if packet.is_up:
                hot_up[slice_index] += electrons_per_packet
            else:
                hot_dn[slice_index] += electrons_per_packet

        hot_tot = [0] * self.num_slices
        mu0_hot_up = [0.0] * self.num_slices
        mu0_hot_dn = [0.0] * self.num_slices

        for i in range(self.num_slices):
            hot_tot[i] = hot_up[i] + hot_dn[i]
            mu0_hot_up[i] = hot_up[i] / (self.slice_length * self.slice_property_list[i].ds_up)
            mu0_hot_dn[i] = hot_dn[i] / (self.slice_length * self.slice_property_list[i].ds_dn)

        mu0_up = [0.0] * self.num_slices
        mu0_dn = [0.0] * self.num_slices

        for i in range(self.num_slices):
            ds_tot_i = self.slice_property_list[i].ds_up + self.slice_property_list[i].ds_dn

            mu0_up[i] = (
                            +self.gamma_list[i] * self.slice_property_list[i].ds_dn -
                            hot_tot[i] / self.slice_length
                        ) / ds_tot_i
            mu0_dn[i] = (
                            -self.gamma_list[i] * self.slice_property_list[i].ds_up -
                            hot_tot[i] / self.slice_length
                        ) / ds_tot_i

        return (self.gamma_list.copy(), mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up, mu0_hot_dn, self.j_hot, self.j_up, self.j_dn)

    def make_ticks(self):
        ticks_slices = []
        for i in range(self.num_slices):
            ticks_slices.append((i + 0.5) * self.slice_length)

        ticks_planes = []
        for i in range(self.num_slices - 1):
            ticks_planes.append((i + 1.0) * self.slice_length)

        return (ticks_slices, ticks_planes)

    def plot(self):
        (ticks_slices, ticks_planes) = self.make_ticks()

        (gamma, mu0_up, mu0_dn, hot_up, hot_dn, mu0_hot_up, mu0_hot_dn, j_hot, j_up, j_dn) = states[n]

        plt.figure(figsize=(9, 9))

        plt.subplot(311)
        plt.plot(ticks_slices, self.gamma_list)
        plt.ylabel("gamma (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.slice_length, -0.1, 0.1])

        plt.subplot(312)
        plt.plot(ticks_slices, mu0_up, ticks_slices, mu0_hot_up)
        plt.ylabel("mu_0 up (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.slice_length, -0.1, 0.1])

        plt.subplot(313)
        plt.plot(ticks_slices, mu0_dn, ticks_slices, mu0_hot_dn)
        plt.ylabel("mu_0 dn (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.slice_length, -0.1, 0.1])

        plt.show()

def make_system():
    system = System()

    system.num_slices = 60
    system.slice_length = 0.5  # (nm)
    system.gamma_list = [0.0] * system.num_slices

    slice_properties_1 = SliceProperties()
    slice_properties_1.ds_up = 70.0  # (eV^-1 nm^-3)
    slice_properties_1.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_1.tau = 10  # (NOT fs)

    slice_properties_2 = SliceProperties()
    slice_properties_2.ds_up = 30.0  # (eV^-1 nm^-3)
    slice_properties_2.ds_dn = 30.0  # (eV^-1 nm^-3)
    slice_properties_2.tau = 100  # (NOT fs)

    plane_properties_1 = PlaneProperties()
    plane_properties_1.alpha_up = 8.0  # (eV^-1 nm^-1 fs^-1)
    plane_properties_1.alpha_dn = 4.0  # (eV^-1 nm^-1 fs^-1)

    plane_properties_2 = PlaneProperties()
    plane_properties_2.alpha_up = 0.4  # (eV^-1 nm^-1 fs^-1)
    plane_properties_2.alpha_dn = 0.4  # (eV^-1 nm^-1 fs^-1)

    system.slice_property_list = [slice_properties_1] * 30
    system.slice_property_list.extend([slice_properties_2] * 30)

    system.plane_property_list = [plane_properties_1] * 30
    system.plane_property_list.extend([plane_properties_2] * 29)

    system.dt = 0.2

    return system
