from .particle import Particle
from .triton import Triton
from scipy import constants as C



# Nucleus of deuterium
class Deuteron(Particle):
    def __init__(self, x=0, y=0, z=0, vx=0, vy=0, vz=0, temperature=298.15,**other_properties):
        Particle.__init__(self, x, y, z, vx, vy, vz, temperature, **other_properties)
        self.mass_MeV = C.value("deuteron mass energy equivalent in MeV")  # [MeV/c2]
        self.mass_kg = C.value("deuteron mass")  # [kg]
        self.charge = C.value("elementary charge")  # [C]
        self.radius = C.value("deuteron rms charge radius")  # [m]
        self.magneticDipoleMoment = C.value("deuteron mag. mom.")  # [J*T**-1]

    def get_color(self):
        return '#ff1a00' # red