from scipy import constants as C


class Laser:
    def __init__(self, wavelength=3.51e-8, intensity=5e18,
                 activation_time=1):  # default values from NIF inertial confinement fusion laser
        self.wavelength = wavelength  # in meters
        self.intensity = intensity  # in W/m^2
        self.frequency = C.speed_of_light / wavelength  # in Hz
        self.activation_time = activation_time

    def get_initial_energy(self, area):
        return self.get_energy(area, self.activation_time)

    def get_energy(self, area, dt):
        photon_energy = C.h * self.frequency  # eV # wavelength in [um]
        n_photons = photon_energy * self.intensity / area * dt
        return C.eV * photon_energy * n_photons * 7.242971666667e+22  # [K]
