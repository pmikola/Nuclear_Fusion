from Particles.helion import Helion
from Particles.triton import Triton
from Particles.neutron import Neutron
from Particles.deuteron import Deuteron
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from chamber import Chamber
from laser import Laser
from matplotlib.widgets import Slider
import random
from matplotlib.offsetbox import AnchoredText
from matplotlib.animation import PillowWriter
import imageio

random.seed(0)
laser = Laser()
matplotlib.use('Qt5Agg')

print('''\nSimulation scenarios:
(1) Chamber full of particles
(2) Pair of particles moving towards each other (with collision)
(3) Pair of particles moving towards each other (without collision)
(4) Pair of particles at an angle (with collision)
(5) Pair of particles at an angle (without collision)
(6) Custom simulation
''')
scenario = int(input("Choose scenario: ") or "1")

chamber = Chamber(laser, scenario=scenario)
stop_time = float(input("Stop Time: ") or "1e-7")

if scenario == 6:
    particle_pairs = int(input("\nNumber of D+T pairs [10]: ") or "10")
    wavelength = float(input("Laser wavelength [3.51e-8 m]: ") or "3.51e-8")
    intensity = float(input("Laser intensity [5e17 W/m^2]: ") or "5e17")
    activation_time = int(input("Laser activation time [1 s]: ") or "1")
    laser = Laser(wavelength, intensity, activation_time)
    chamber = Chamber(laser, scenario=scenario, particle_pairs=particle_pairs)

chamber.create_particles()

fig = plt.figure()

grid = plt.GridSpec(4, 4, wspace=6, hspace=0.6)
at = fig.add_subplot(grid[:2, :2])
ax = fig.add_subplot(grid[2:, :2], projection='3d')
aa = fig.add_subplot(grid[0, 2:])
ay = fig.add_subplot(grid[1, 2:])
az = fig.add_subplot(grid[2, 2:])
ab = fig.add_subplot(grid[3, 2:])

ax.set_xlim([0, chamber.x])
ax.set_ylim([0, chamber.y])
ax.set_zlim([0, chamber.z])
ax.set_xlabel('X axis [m] ', fontsize=8, labelpad=4)
ax.set_ylabel('Y axis [m] ', fontsize=8, labelpad=4)
ax.set_zlabel('Z axis [m] ', fontsize=8, labelpad=4)
ax.xaxis.set_tick_params(labelsize=6)
ax.yaxis.set_tick_params(labelsize=6)
ax.zaxis.set_tick_params(labelsize=6)

boax = ax.get_position()
boax.y0 = boax.y0 + 0.1
boax.y1 = boax.y1 + 0.1
boax.x0 = boax.x0 + 0
boax.x1 = boax.x1 + 0
ax.set_position(boax)

ay.set_ylabel('Temperature [ K ] ', fontsize=8)
ay.grid(True)
# ay.yaxis.set_ticks_position("right")
ay.yaxis.set_label_position("right")
# ay.yaxis.tick_right()


az.set_ylabel('Pressure [ Pa ] ', fontsize=8)
az.grid(True)
az.yaxis.set_label_position("right")

aa.set_ylabel('Sub Energy [ J ] ', fontsize=8)
aa.grid(True)
aa.yaxis.set_label_position("right")

ab.set_xlabel('Time [ s ]')
ab.set_ylabel('Avg Velocity [m/s]', fontsize=8)
ab.grid(True)
ab.yaxis.set_label_position("right")

at.get_yaxis().set_visible(False)
at.get_xaxis().set_visible(False)
at.axis("off")
box = at.get_position()
box.y0 = box.y0 + 0.00
box.y1 = box.y1 + 0.00
box.x0 = box.x0 + 0.06
box.x1 = box.x1 + 0.06
at.set_position(box)

lines, x, y, z, helionx, heliony, helionz, \
deuteronx, deuterony, deuteronz, tritonx, tritony, tritonz, \
neutronx, neutrony, neutronz, Time, Velocity, \
Temperature, SubEnergy, ReactionCount, TotalEnergy, Pressure, title = ([] for i in range(24))
helion_color = Helion().get_color()
neutron_color = Neutron().get_color()
triton_color = Triton().get_color()
deuteron_color = Deuteron().get_color()
laserEnergy = chamber.laser_energy
Dt = chamber.dt

while True:
    chamber.update_particles()
    for p in chamber.particles:
        if isinstance(p, Helion):
            helionx.append(p.position.x)
            heliony.append(p.position.y)
            helionz.append(p.position.z)
        elif isinstance(p, Deuteron):
            deuteronx.append(p.position.x)
            deuterony.append(p.position.y)
            deuteronz.append(p.position.z)
        elif isinstance(p, Triton):
            tritonx.append(p.position.x)
            tritony.append(p.position.y)
            tritonz.append(p.position.z)
        else:
            neutronx.append(p.position.x)
            neutrony.append(p.position.y)
            neutronz.append(p.position.z)

    Time.append(chamber.Time)
    SubEnergy.append(chamber.sub_energy)
    ReactionCount.append(chamber.reaction_count)
    TotalEnergy.append(chamber.total_energy_released)
    Temperature.append(chamber.Temperature)
    Pressure.append(chamber.Pressure)
    Velocity.append(chamber.avg_vel)

    lineh, = ax.plot(helionx, heliony, helionz, color='None', marker='.', markeredgecolor=helion_color,
                     markerfacecolor=helion_color)
    lined, = ax.plot(deuteronx, deuterony, deuteronz, color='None', marker='.', markeredgecolor=deuteron_color,
                     markerfacecolor=deuteron_color)
    linet, = ax.plot(tritonx, tritony, tritonz, color='None', marker='.', markeredgecolor=triton_color,
                     markerfacecolor=triton_color)
    linen, = ax.plot(neutronx, neutrony, neutronz, color='None', marker='.', markeredgecolor=neutron_color,
                     markerfacecolor=neutron_color)
    line2, = ay.semilogy(Time, Temperature, color='red', marker=',')
    line3, = az.semilogy(Time, Pressure, color='green', marker=',')
    line4, = aa.semilogy(Time, SubEnergy, color='blue', marker=',')
    line5, = ab.semilogy(Time, Velocity, color='purple', marker=',')

    anchored_text = AnchoredText(
        "Laser Init Energy:           " + '{:<.3e}'.format(laserEnergy * 1.38064878066852e-23 / Dt) + "[W]" + "\n"
        + "Particle Number:            " + str(chamber.particle_pairs * 2) + "\n"
        + "Volume:                         " + '{:<.2e}'.format(chamber.x * chamber.y * chamber.z) + "[m3]" + "\n"
        + "Time:                             " + '{:<.3e}'.format(chamber.Time) + " [s]" + "\n"
        + "Subtracted Energy:        " + '{:<.3e}'.format(chamber.sub_energy) + " [J] " + "\n"
        + "Total reactions:              " + str(chamber.reaction_count) + "\n"
        + "Total released energy:   " + '{:<.2f}'.format(
            chamber.total_energy_released) + " [MeV] " + "\n"
        + "Temperature:                 " + '{:<.3e}'.format(
            chamber.Temperature) + " [K]" + "\n"
        + "Pressure:                       " + '{:<.3e}'.format(chamber.Pressure) + " [Pa]",
        loc=5, prop=dict(size=8))

    title = at.add_artist(anchored_text)

    lines.append([lineh, lined, linet, linen, title, line2, line3, line4, line5])
    x.clear(), y.clear(), z.clear(), helionx.clear(), heliony.clear(), helionz.clear(), deuteronx.clear(), deuterony.clear(), deuteronz.clear()
    tritonx.clear(), tritony.clear(), tritonz.clear(), neutronx.clear(), neutrony.clear(), neutronz.clear()
    # if chamber.Time <= stop_time:
    if chamber.Time <= stop_time:
        continue
    else:
        break

# file_name = "scenario" + str(scenario)
# file_name = "./" + file_name + '.gif'
ani = animation.ArtistAnimation(fig, lines, interval=10, blit=True)
# ani.save(file_name, writer='pillow', fps=30)
figManager = plt.get_current_fig_manager()
figManager.window.resize(700, 700)

plt.show()
