from math import sqrt
from scipy import constants as C
from vectormath import Vector3

max_speed = 1e8

class Particle:
    def __init__(self, x=0, y=0, z=0, vx=0, vy=0, vz=0, temperature=298.15, **other_properties):
        self.position = Vector3(x, y, z)
        self.velocity = Vector3(vx, vy, vz)  # [m/s]
        self.temperature = temperature  # [K]

        for name, value in other_properties.items():
            setattr(self, name, value)

    def get_speed(self):
        v = self.velocity
        return sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2)

    # modify particle velocity based on influence from some particle
    def get_influence_from(self, particle, dt):
        dist = self.get_distance_to(particle)
        if dist > 0:
            force = (C.k * self.charge * particle.charge) / (dist ** 2)
            acc = force / self.mass_kg
            additional_velocity = dt * acc * -self.get_vector_to(particle)
            self.velocity += additional_velocity
            if self.velocity.length > max_speed:
                self.velocity = self.velocity.normalize() * max_speed

    def get_distance_to(self, particle):
        return (particle.position - self.position).length

    # modify position by current velocity value
    def update_position(self, dt):
        vel = self.velocity * dt
        self.position += vel

    def add_energy(self, energy):
        current_energy = self.get_kinetic_energy()
        new_speed = sqrt(2 * (current_energy + energy) / self.mass_kg)
        self.velocity.normalize()
        self.velocity *= new_speed
        
    # get unique color for this particle type (used for plotting)
    def get_color(self):
        return '#000000'

    # get normalized vector from self to particle
    def get_vector_to(self, particle):
        return (particle.position - self.position).normalize()

    def get_kinetic_energy(self):
        return self.mass_kg * self.get_speed() ** 2 / 2
