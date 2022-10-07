from math import sqrt
from random import Random

from Particles import particle
from Particles.particle import Particle
from Particles.triton import Triton
from Particles.deuteron import Deuteron
from Particles.helion import Helion
from Particles.neutron import Neutron
from scipy import constants as C
from scipy.stats import norm
import random

energy_released_in_MeV = 17.59
neutron_energy_ratio = 0.7987
MeV_in_Joules = 1.6021773E-13
start_speeds = [1e6, 5e-12, 2e-12, 5e-12, 2e-12]
wallCoeff = 0.3  # percent of grabbed energy
chamber_sizes = [1e-2, 2e-10, 2e-10, 2e-10, 2e-10]
fusion_distance_ratio = 5e-2


class Chamber:
    def __init__(self, laser, x=1, y=1, z=1, scenario=1, particle_pairs=5):
        # chamber dimensions
        self.laser = laser
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.particle_pairs = particle_pairs
        self.scenario = scenario
        self.particles = list()
        self.total_energy_released = 0
        self.reaction_count = 0
        self.sub_energy = 0
        self.Temperature = 298.15  # [ K ]
        self.Pressure = 0  # [ Pa ]
        self.KEavg = 0
        self.Time = 0
        self.iterator = 1
        self.dt = 1e-9  # [ s ]
        self.min_fusion_dist = 5e-4
        self.set_up_scenario()
        self.surface_area = self.get_chamber_surface_area()
        self.laser_energy = self.laser.get_initial_energy(self.surface_area)
        self.resultant_velocity = 0
        self.avg_velocity = 0
        self.avg_vel = 0

    def set_up_scenario(self):
        if self.scenario <= 5:
            chamber_size = chamber_sizes[self.scenario - 1]
            self.x = self.y = self.z = chamber_size
            self.min_fusion_dist = chamber_size * fusion_distance_ratio

    # create particles with random positions and velocities
    def create_particles(self):
        if self.scenario <= 5:
            start_speed = start_speeds[self.scenario - 1]
        else:
            particle_energy = self.laser_energy / (self.particle_pairs * 2)
            start_speed = sqrt(2 * particle_energy / Deuteron().mass_kg)
        if self.scenario == 1 or self.scenario == 6:
            for _ in range(0, self.particle_pairs):
                x, y, z = self.get_random_position()
                vx, vy, vz = self.get_random_velocity(1)
                deuteron = Deuteron(x, y, z, vx, vy, vz)
                deuteron.velocity = deuteron.velocity.normalize() * start_speed
                self.particles.append(deuteron)
                x, y, z = self.get_random_position()
                vx, vy, vz = self.get_random_velocity(0)
                triton = Triton(x, y, z, vx, vy, vz)
                triton.velocity = triton.velocity.normalize() * start_speed
                self.particles.append(triton)
        elif self.scenario == 2 or self.scenario == 3:
            deuteron = Deuteron(0, self.y / 2, self.z / 2, start_speed)
            triton = Triton(self.x, self.y / 2, 1.01 * self.z / 2, -start_speed)
            self.particles.append(deuteron)
            self.particles.append(triton)
        elif self.scenario == 4 or self.scenario == 5:
            deuteron = Deuteron(0, 0, self.z / 2, start_speed, start_speed)
            triton = Triton(self.x, 0, 1.01 * self.z / 2, -start_speed, start_speed)
            self.particles.append(deuteron)
            self.particles.append(triton)

    # update particle parameters
    def update_particles(self):
        current_max_v = 0
        for i in range(0, len(self.particles)):
            p1 = self.particles[i]
            for j in range(0, len(self.particles)):
                if i != j:
                    p2 = self.particles[j]
                    if self.fusion_can_occur(p1, p2):
                        helion, neutron = self.execute_fusion(p1, p2)
                        self.particles[i] = helion
                        self.particles[j] = neutron
                    p1.get_influence_from(p2, self.dt)
            speed = p1.get_speed()
            if speed > current_max_v:
                current_max_v = speed
        if self.scenario > 0 and current_max_v != 0:
            self.set_dt(current_max_v)
        for p in self.particles:
            p.update_position(self.dt)

        self.clip_to_bounds()
        self.Time += self.dt

    def fusion_can_occur(self, p1, p2):
        if p1.get_distance_to(p2) < self.min_fusion_dist:
            if (isinstance(p1, Triton) and isinstance(p2, Deuteron)) or (
                    isinstance(p1, Deuteron) and isinstance(p2, Triton)):
                return True
        return False

    # carry out nuclear fusion
    def execute_fusion(self, deuteron, triton):
        self.reaction_count += 1
        self.total_energy_released += energy_released_in_MeV

        helion = Helion(triton.position.x, triton.position.y, triton.position.z)
        neutron = Neutron(deuteron.position.x, deuteron.position.y, deuteron.position.z)

        input_energy = triton.get_kinetic_energy() + deuteron.get_kinetic_energy()
        output_energy = input_energy + energy_released_in_MeV * MeV_in_Joules
        neutron_energy = neutron_energy_ratio * output_energy  # in Joules
        helion_energy = output_energy - neutron_energy  # in Joules

        helion_speed = sqrt(2 * helion_energy / helion.mass_kg)
        helion.velocity = triton.velocity.normalize() * helion_speed

        neutron_speed = sqrt(2 * neutron_energy / neutron.mass_kg)
        neutron.velocity = deuteron.velocity.normalize() * neutron_speed

        return helion, neutron

    # clip particles inside the chamber
    def clip_to_bounds(self):
        i = 0
        self.KEavg = 0
        for particle in self.particles:
            m = particle.mass_kg
            p = particle.position
            v = particle.velocity
            if p.x < 0:
                v.x = self.set_bound_damping(v.x, m)
                v.x = -v.x
                p.x = 0
            elif p.x > self.x:
                v.x = self.set_bound_damping(v.x, m)
                v.x = -v.x
                p.x = self.x
            if p.y < 0:
                v.y = self.set_bound_damping(v.y, m)
                v.y = -v.y
                p.y = 0
            elif p.y > self.y:
                v.y = self.set_bound_damping(v.y, m)
                v.y = -v.y
                p.y = self.y
            if p.z < 0:
                v.z = self.set_bound_damping(v.z, m)
                v.z = -v.z
                p.z = 0
            elif p.z > self.z:
                v.z = self.set_bound_damping(v.z, m)
                v.z = -v.z
                p.z = self.z
            i += 1
            self.resultant_velocity = sqrt(v.x ** 2 + v.y ** 2 + v.z ** 2)
            self.KEavg += ((m * (self.resultant_velocity ** 2)) / 2)
            self.iterator += 1
            self.avg_velocity += self.resultant_velocity
        self.avg_vel = self.avg_velocity / len(self.particles)
        self.Pressure = (len(self.particles) * self.KEavg / (self.x * self.y * self.z))
        self.Temperature = ((self.KEavg / len(self.particles)) / (3 * C.k))
        self.avg_velocity = 0

    def get_random_position(self):
        random = Random()
        x = random.uniform(0, self.x)
        y = random.uniform(0, self.y)
        z = random.uniform(0, self.z)
        return x, y, z

    def gas_velocity(self, flag):
        vsqr = 0
        if flag == 1:
            m = C.value("deuteron mass")
        if flag == 0:
            m = C.value("triton mass")
        T = self.Temperature
        vsqr += ((3 / 2) * C.k * T) / (1 * m / 2)
        return sqrt(vsqr) / self.particle_pairs * 2

    def get_random_velocity(self, flag):
        random = Random()
        v = self.gas_velocity(flag)
        vx = random.uniform(-v, v)
        vy = random.uniform(-v, v)
        vz = random.uniform(-v, v)
        return vx, vy, vz

    def set_bound_damping(self, v, m):
        KELost = (m * v ** 2 / 2) * wallCoeff
        self.sub_energy += KELost
        particle_count = len(self.particles)
        energy_per_particle = KELost / particle_count
        for i in range(0, particle_count):
            p = self.particles[i]
            p.add_energy(energy_per_particle)
        return v * sqrt(1 - wallCoeff)

    def set_dt(self, max_v):
        if self.scenario == 1 or self.scenario == 6:
            n_frames = 5
        else:
            n_frames = 20
        self.dt = (self.x / n_frames) / max_v  # so that it takes n frames to move across chamber

    def get_chamber_surface_area(self):
        return 2 * (self.x * self.y + self.x * self.z + self.y * self.z)
