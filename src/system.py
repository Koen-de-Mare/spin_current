import matplotlib.pyplot as plt


from src.material import Material


class System:
    def __init__(self):
        # state
        self.gamma_list = []
        self.t = 0.0

        # properties
        self.num_slices = 0
        self.slice_length = 0.0  # (nm)
        self.material = Material()
        self.dt = 0.0  # (fs)

    def step(self):
        self.t += self.dt
        xi = self.material.xi()

        # "current"
        j_list = [0.0] * (self.num_slices - 1)

        for i in range(self.num_slices - 1):
            j_list[i] = - xi * (self.gamma_list[i+1] - self.gamma_list[i]) / self.slice_length;

        # time derivative of gamma
        gamma_prime = [0.0] * self.num_slices
        for i in range(1, self.num_slices):
            gamma_prime[i] += j_list[i - 1] / self.slice_length
        for i in range(0, self.num_slices - 1):
            gamma_prime[i] -= j_list[i] / self.slice_length
        for i in range(0, self.num_slices):
            gamma_prime[i] -= \
                self.gamma_list[i] * (1.0 / self.material.ds_up + 1.0 / self.material.ds_dn) / self.material.tau

        for i in range(0, self.num_slices):
            self.gamma_list[i] += self.dt * gamma_prime[i]

    def make_ticks(self):
        ticks = []
        for i in range(self.num_slices):
            ticks.append(i * self.slice_length)
        return ticks

    def plot(self):
        ticks = self.make_ticks()
        plt.plot(ticks, self.gamma_list)
        plt.xlabel("depth (nm)")
        plt.ylabel("gamma (eV)")
        plt.show()

    def stability(self) -> float:
        # this number should be smaller than 1 and larger than 0 for hope of accurate results
        return 1 - 2 * self.dt * self.material.xi() / (self.slice_length * self.slice_length)


def make_system():
    system = System()

    system.num_slices = 100
    system.slice_length = 1.0  # (nm)
    system.gamma_list = [0.0] * system.num_slices

    system.material.ds_up = 90.0  # (eV^-1 nm^-3)
    system.material.ds_dn = 10.0  # (eV^-1 nm^-3)
    system.material.alpha_up = 8.0  # (eV^-1 nm^-1 fs^-1)
    system.material.alpha_dn = 2.0  # (eV^-1 nm^-1 fs^-1)
    system.material.tau = 100000  # (NOT fs)

    system.dt = 1.0

    return system
