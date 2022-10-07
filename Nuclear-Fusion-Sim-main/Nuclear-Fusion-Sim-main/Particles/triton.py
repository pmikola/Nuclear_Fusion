from .particle import Particle
from scipy import constants as C

# Nucleus of tritium
class Triton(Particle):
    def __init__(self, x=0, y=0, z=0, vx=0, vy=0, vz=0, temperature=298.15,**other_properties):
        Particle.__init__(self, x, y, z, vx, vy, vz, temperature, **other_properties)
        self.mass_MeV = C.value("triton mass energy equivalent in MeV")  # [MeV/c2]
        self.mass_kg = C.value("triton mass")  # [kg]
        self.charge = 2 * C.value("elementary charge")  # [C]
        self.magneticDipoleMoment = C.value("triton mag. mom.")  # [J*T**-1]

    def get_color(self):
        return '#2254dd' # blue