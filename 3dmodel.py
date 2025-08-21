import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque
import time

matplotlib.use("Agg")

"""Globals"""

# constants
c = 299792458  # m/s
n = 1.58       # EJ-200 refractive index
v = c / n

# empty sets
fading_particles = []
dead_particles = []

"""Functions"""

class Particle:
    def __init__(self, position, velocity, axes):
        # position/velocity are 3D now: [x, y, z]
        self.position = np.array(position, dtype=float)
        self.velocity = np.array(velocity, dtype=float)
        self.velocity *= v / np.linalg.norm(self.velocity) #moved this up. Velocity does not exceed max value

        ax_xy, ax_xz, ax_yz = axes
        # dots for each projection
        self.dot_xy, = ax_xy.plot([], [], 'o', markersize=4, color='blue')
        self.dot_xz, = ax_xz.plot([], [], 'o', markersize=4, color='blue')
        self.dot_yz, = ax_yz.plot([], [], 'o', markersize=4, color='blue')

        self.trail_length = trail_length
        self.positions = deque(maxlen=trail_length)
        # trails per projection
        self.trail_xy, = ax_xy.plot([], [], '-', lw=1, color='blue')
        self.trail_xz, = ax_xz.plot([], [], '-', lw=1, color='blue')
        self.trail_yz, = ax_yz.plot([], [], '-', lw=1, color='blue')

        self.distance_traveled = 0.0  # in cm

        self.dead_ring_size = 4
        self.dead_ring_fade_rate = 0.125

        self.absorbed = False
        self.counted = False

    def pmt_check(self):
        for p in pmts:
            if (p['xmin'] <= self.position[0] <= p['xmax'] and #now this if statement asks if inside the volume of the PMT
                p['ymin'] <= self.position[1] <= p['ymax'] and
                p['zmin'] <= self.position[2] <= p['zmax']):
                self.counted = True
                pmt_hits[p['channel'] - 1] += 1 #add a count to the respective PMT

                sumtop = pmt_hits[0] + pmt_hits[1]
                sumbot = pmt_hits[2] + pmt_hits[3]
                print(f"Top sum hits: {sumtop}")
                print(f"Bot sum hits: {sumbot}")

                break

    def update_active(self, dt, bounds): #update every frame
        if self.absorbed:
            return

        dx = self.velocity * dt * 100.0  # x100 to make m/s -> cm/s

        # attenuation
        self_absorption_prob = np.exp(-np.linalg.norm(dx) / attenuation)
        if np.random.rand() > self_absorption_prob:
            self.absorbed = True
            return

        # step forward
        self.position += dx

        # walls in 3D
        for i in range(3):  # 0=x, 1=y, 2=z
            # die if hit right wall
            if i == 0 and self.position[0] >= bounds[0]:
                self.absorbed = True
                return

            if not (0 < self.position[i] < bounds[i]): #if outside slab
                normal = np.zeros(3) #set up a normal vector
                normal[i] = -np.sign(self.velocity[i]) #assigns inward pointing normal
                perturbed_normal = normal + np.random.normal(0,max_pert,3)
                perturbed_normal /= np.linalg.norm(perturbed_normal) #renormalize

                # incidence angle relative to normal of wall i
                v_normal = np.dot(self.velocity,perturbed_normal) #component of photon velocity parallel to normal
                v_parallel = np.linalg.norm(np.cross(self.velocity, perturbed_normal)) #component of photon velocity parallel to wall
                angle = np.arctan2(v_parallel,abs(v_normal))

                if angle < TIR:
                    self.pmt_check()
                    self.absorbed = True
                    return
                else:
                    # reflect and renormalize speed
                    self.position[i] = np.clip(self.position[i], 0, bounds[i]) #put us back inside the slab so we dont get stuck outside
                    self.velocity = self.velocity - 2*v_normal*perturbed_normal #reflect with pertubation
                    self.velocity *= v / np.linalg.norm(self.velocity) #make sure we dont exceed our maximum speed


            self.positions.append(self.position.copy())
            x, y, z = self.position
            xs, ys, zs = np.array(self.positions).T


        self.trail_xy.set_data(xs, ys) #3 planes -> 3 plots -> 3 trails
        self.trail_xz.set_data(xs, zs)
        self.trail_yz.set_data(ys, zs)

        self.dot_xy.set_data([x], [y]) #same for dots
        self.dot_xz.set_data([x], [z])
        self.dot_yz.set_data([y], [z])

    def update_fading(self):
        edge_color = 'green' if self.counted else 'red'
        if self.dead_ring_size > 0:
            self.dead_ring_size = max(0, self.dead_ring_size - self.dead_ring_fade_rate)
            for d in (self.dot_xy, self.dot_xz, self.dot_yz):
                d.set_markersize(self.dead_ring_size)
                d.set_markeredgecolor(edge_color)
                d.set_markerfacecolor('none')

            if self.positions:
                self.positions.popleft()

            if self.positions:
                xs, ys, zs = np.array(self.positions).T
            else:
                xs = ys = zs = np.array([])

            self.trail_xy.set_data(xs, ys)
            self.trail_xz.set_data(xs, zs)
            self.trail_yz.set_data(ys, zs)
        return

# animation functions
def init():
    artists = []
    for p in active_particles + fading_particles:
        # XY
        p.dot_xy.set_data([p.position[0]], [p.position[1]])
        p.trail_xy.set_data([], [])
        # XZ
        p.dot_xz.set_data([p.position[0]], [p.position[2]])
        p.trail_xz.set_data([], [])
        # YZ
        p.dot_yz.set_data([p.position[1]], [p.position[2]])
        p.trail_yz.set_data([], [])
        artists.extend([p.dot_xy, p.trail_xy, p.dot_xz, p.trail_xz, p.dot_yz, p.trail_yz])
    return artists

def animate(frame):
    global lastUpdate, active_particles, fading_particles, dead_particles

    dt = 5.3e-12 * 5 #this dt makes each frame approximately a 0.5 cm step

    
    for p in active_particles[:]:
        p.update_active(dt, dims)
        if p.absorbed:
            active_particles.remove(p)
            fading_particles.append(p)

    for p in fading_particles[:]:
        p.update_fading()
        if p.dead_ring_size == 0:
            fading_particles.remove(p)
            dead_particles.append(p)

    sim_time = frame * dt  # seconds
    timer_text.set_text(f"Time: {sim_time*1e9:.2f} ns") #not sure if this is actually accurate
    alive_text.set_text(f"Alive: {len(active_particles)}")
    pmt1_text.set_text(f"CH1: {pmt_hits[0]}")
    pmt2_text.set_text(f"CH2: {pmt_hits[1]}")
    pmt3_text.set_text(f"CH3: {pmt_hits[2]}")
    pmt4_text.set_text(f"CH4: {pmt_hits[3]}")
    sumtop_text.set_text(f"Top sum: {pmt_hits[0] + pmt_hits[3]}")
    sumbot_text.set_text(f"Bot sum: {pmt_hits[1] + pmt_hits[3]}")

    artists = []
    for p in active_particles + fading_particles:
        artists.extend([p.dot_xy, p.trail_xy, p.dot_xz, p.trail_xz, p.dot_yz, p.trail_yz])
    return artists

def draw_pmts(ax_xy, ax_xz, ax_yz, color_map):
    # XY: just draw the squares as before
    for p in pmts:
        width = p['xmax'] - p['xmin']
        height = p['ymax'] - p['ymin']
        ax_xy.add_patch(plt.Rectangle((p['xmin'], p['ymin']), width, height,color=color_map[p['channel']], alpha=0.6))
        ax_xy.text((p['xmin']+p['xmax'])/2, (p['ymin']+p['ymax'])/2, str(p['channel']),fontsize=10, fontweight='bold', ha='center', va='center')

    # XZ: pmts overlap in this plane, so it looks purple. might just remove this drawing entirely
    for p in pmts:
        width = p['xmax'] - p['xmin']
        depth = p['zmax'] - p['zmin']
        ax_xz.add_patch(plt.Rectangle((p['xmin'], p['zmin']), width, depth,color=color_map[p['channel']], alpha=0.3))

    # YZ: same color pmts overlap now, but they show up cleanly as boxes on either side of the slab
    for p in pmts:
        height = p['ymax'] - p['ymin']
        depth = p['zmax'] - p['zmin']
        ax_yz.add_patch(plt.Rectangle((p['ymin'], p['zmin']), height, depth,color=color_map[p['channel']], alpha=0.3))

def plotter(dims):
    # three plots
    fig, (ax_xy, ax_xz, ax_yz) = plt.subplots(3, 1, figsize=(8,12))
    fig.subplots_adjust(left=0.15, hspace=0.25)

    # XY box
    ax_xy.set_xlim(-10, dims[0] + 10)
    ax_xy.set_ylim(-10, dims[1] + 20)
    ax_xy.set_aspect('equal')
    ax_xy.set_title('XY')
    ax_xy.add_patch(plt.Rectangle((0, 0), dims[0], dims[1], fill=False))

    # XZ box
    ax_xz.set_xlim(-10, dims[0] + 10)
    ax_xz.set_ylim(-5, dims[2] + 5)
    ax_xz.set_aspect('equal')
    ax_xz.set_title('XZ')
    ax_xz.add_patch(plt.Rectangle((0, 0), dims[0], dims[2], fill=False))

    # YZ box
    ax_yz.set_xlim(-10, dims[1] + 20)
    ax_yz.set_ylim(-5, dims[2] + 5)
    ax_yz.set_aspect('equal')
    ax_yz.set_title('YZ')
    ax_yz.add_patch(plt.Rectangle((0, 0), dims[1], dims[2], fill=False))

    # Color map: top (CH1, CH2) red; bottom (CH3, CH4) blue
    color_map = {1: 'red',
                 #2: 'blue', #this is here for the 2 PMT case (broida geometry) just uncomment this and comment the others out 
                 2: 'red', 
                 3: 'blue', 
                 4: 'blue'
                 }
    draw_pmts(ax_xy, ax_xz, ax_yz, color_map)

    # HUD

    timer_text = ax_xy.text(0.02, 1.02, "", transform=ax_xy.transAxes, fontsize=10, va="bottom", ha="left")
    alive_text = ax_xy.text(0.02, 0.96, "", transform=ax_xy.transAxes, fontsize=10, va="top", ha="left")

    pmt1_text = ax_xy.text(0.38, 1.02, "", transform=ax_xy.transAxes, fontsize=10, va="bottom", ha="left", color="red")
    pmt2_text = ax_xy.text(0.52, 1.02, "", transform=ax_xy.transAxes, fontsize=10, va="bottom", ha="left", color="red")
    pmt3_text = ax_xy.text(0.66, 1.02, "", transform=ax_xy.transAxes, fontsize=10, va="bottom", ha="left", color="blue")
    pmt4_text = ax_xy.text(0.80, 1.02, "", transform=ax_xy.transAxes, fontsize=10, va="bottom", ha="left", color="blue")
    sumtop_text = ax_xy.text(0.38, 0.96, "", transform=ax_xy.transAxes, fontsize=10, va="top", ha="left", color="red")
    sumbot_text = ax_xy.text(0.66, 0.96, "", transform=ax_xy.transAxes, fontsize=10, va="top", ha="left", color="blue")

    return fig, (ax_xy, ax_xz, ax_yz), timer_text, alive_text, pmt1_text, pmt2_text, pmt3_text, pmt4_text, sumtop_text, sumbot_text

#randomized function
def particlegen(n, pos, axes):
    particles = []
    for _ in range(n):
        # random direction on a sphere
        phi = np.random.uniform(0, 2*np.pi) #azimuthal angle
        costheta = np.random.uniform(-1, 1) #cos of the polar angle
        sintheta = np.sqrt(1 - costheta**2) #sin of the polar angle using sin^2 + cos^2 = 1
        direction = np.array([np.cos(phi) * sintheta, np.sin(phi) * sintheta, costheta])
        velocity = v * direction
        particles.append(Particle(pos, velocity, axes))
    return particles

"""main"""

#objects
dims = np.array([110.0, 76.0, 5.0])
pmts = [
    {'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': 76,  'ymax': 81,  'zmin': 0, 'zmax': 5, 'channel': 1},  
    {'xmin': 60-2.5, 'xmax': 60+2.5, 'ymin': 76,  'ymax': 81,  'zmin': 0, 'zmax': 5, 'channel': 3}, 
    {'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': -5,  'ymax': 0,   'zmin': 0, 'zmax': 5, 'channel': 2}, 
    {'xmin': 60-2.5, 'xmax': 60+2.5, 'ymin': -5,  'ymax': 0,   'zmin': 0, 'zmax': 5, 'channel': 4}, 
]

#physics
TIR = np.radians(39)  
attenuation = 380.0 # cm
max_pert = 0.03 #standard deviation of surface from being perfectly flat

#plotting and tracking
trail_length = 1 #doesnt really show up this way. may not be worth having for more than like 10 particles
pmt_hits = [0] * len(pmts)

# 3 plane figure
fig, axes, timer_text, alive_text, pmt1_text, pmt2_text, pmt3_text, pmt4_text, sumtop_text, sumbot_text = plotter(dims) #use this for 4 pmts
#fig, axes, timer_text, alive_text, pmt1_text, pmt2_text = plotter(dims) #use this for 2 pmts 

# make the particles
spawn = [35, 40, dims[2]/2]
active_particles = particlegen(10, spawn, axes)

start_time = time.time()
lastUpdate = time.time()
animation = FuncAnimation(fig, animate, interval=10, frames=500, blit=True, init_func=init)
animation.save("particle_sim3d2.gif", writer='pillow', fps=10)
end_time = time.time()

print(f"{end_time - start_time:.2f}")
