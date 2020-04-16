class Material:
    def __init__(self):
        self.ds_up = 0
        self.ds_dn = 0
        self.alpha_up = 0
        self.alpha_dn = 0
        self.tau = 0

    def xi(self) -> float:
        return self.alpha_up * self.alpha_dn * (self.ds_up / self.ds_dn + self.ds_dn / self.ds_up + 2.0) / ((self.ds_up + self.ds_dn) * (self.alpha_up + self.alpha_dn))
