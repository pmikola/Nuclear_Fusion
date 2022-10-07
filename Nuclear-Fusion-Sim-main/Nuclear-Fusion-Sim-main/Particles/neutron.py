from .particle import Particle
from scipy import constants as C

class Neutron(Particle):
    def __init__(self, x=0, y=0, z=0, vx=0, vy=0, vz=0, temperature=298.15,**other_properties):
        Particle.__init__(self, x, y, z, vx, vy, vz, temperature, **other_properties)
        self.type = "Baryon"
        self.statistics = "Fermion"
        self.compounds = {"u", "d", "d"}
        self.mass_MeV = C.value("neutron mass energy equivalent in MeV")  # [MeV/c2]
        self.mass_kg = C.value("neutron mass")  # [kg]
        self.charge = 0  # [C]
        self.spin = 1 / 2
        self.radius = 0.8 * 10 ** -15  # [m]
        self.meanLifetime = 879.466  # [s]
        self.magneticDipoleMoment = C.value("neutron mag. mom.")  # [J*T**-1]
        self.electricDipoleMoment = 2.9 * 10 ** -26  # [e*cm]
        self.electricPolarizability = 1.1615 * 10 ** -3  # [f*m**3]
        self.magneticPolarizability = 3.720 * 10 ** -4  # [f*m**3]

    def get_color(self):
        return '#28d73c' # green