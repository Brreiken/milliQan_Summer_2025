import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque
import time
import os
import csv

matplotlib.use("Agg")

"""Globals"""

#constants
c = 299792458 #speed of light
n = 1.58 #refractive index of ej200 plastic scintillator
v = c/n #our reduced speed for traveling through a material

#empty sets
pmt_hits=[]
fading_particles = []
dead_particles = []

"""Functions"""

class Particle:
    def __init__(self,position,velocity, ax):
        self.position = np.array(position,dtype=float) #basic physics parameters
        self.velocity = np.array(velocity,dtype=float)
        self.dot, = ax.plot([],[],'o', markersize=4,color='blue') #and a graphic representation

        self.trail_length = trail_length
        self.positions = deque(maxlen=trail_length) #an array of previous positions that will only be as long as trail length
        self.trail_line, = ax.plot([],[],'-',lw=1,color='blue') 
        self.distance_traveled = 0.0  

        self.dead_ring_size = 4  
        self.dead_ring_fade_rate = 0.125  

        self.absorbed = False 
        self.counted = False

    def pmt_check(self):
        for p in pmts:
            if (p['xmin'] <= self.position[0] <= p['xmax'] and p['ymin'] <= self.position[1] <= p['ymax']): #if in PMT
                self.counted = True
                pmt_hits[p['channel'] - 1] += 1 #count corresponding to correct channel
                break

    def update_active(self,dt,bounds): 
        if self.absorbed:
            return

        dx = self.velocity*dt*100 
        self.position += dx

        self_absorption_prob = np.exp(-np.linalg.norm(dx)/attenuation) #attenuation/self absorption check
        if np.random.rand() > self_absorption_prob:
            self.absorbed=True
            return

        for i in range(2):

            if i == 0 and self.position[0] >= bounds[0]: #kill wall, comment out if you dont want it
                self.absorbed = True
                return

            if not 0 < self.position[i] < bounds[i]:
                self.velocity[i] *= 1 + np.random.uniform(-0.04,0.04) #perturbate the velocity before we check for total internal reflection
                self.velocity[i-1] *= 1 + np.random.uniform(-0.04,0.04) #it effectively simulates surface roughness
                angle = np.arctan2(abs(self.velocity[1-i]),abs(self.velocity[i]))
                if angle < TIR: 
                    self.pmt_check() #if absorbed into wall, check whether its in the PMT or not
                    self.absorbed = True
                    return
                else:
                    self.position[i] = np.clip(self.position[i],0,bounds[i]) #get back in the box
                    self.velocity[i] *= -1 #reflect
                    self.velocity *= v/np.linalg.norm(self.velocity) #renormalize to not exceed speed of light

        self.positions.append(self.position.copy()) 
        x,y = zip(*self.positions)
        self.trail_line.set_data(x,y)
        self.dot.set_data([self.position[0]], [self.position[1]])

    def update_fading(self): #graphic stuff
        edge_color = 'green' if self.counted else 'red'
        if self.dead_ring_size > 0:
            self.dead_ring_size = max(0, self.dead_ring_size - self.dead_ring_fade_rate)
            self.dot.set_markersize(self.dead_ring_size)
            self.dot.set_markeredgecolor(edge_color)
            self.dot.set_markerfacecolor('none')

            if self.positions:
                self.positions.popleft()
            x, y = zip(*self.positions) if self.positions else ([], [])
            self.trail_line.set_data(x, y)
        return

#animation functions

def init():
    artists = [] 
    for p in active_particles + fading_particles:
        p.dot.set_data([p.position[0]], [p.position[1]])
        p.trail_line.set_data([], [])
        artists.extend([p.dot, p.trail_line])
    return artists

def animate(frame): 
    global lastUpdate, active_particles, fading_particles, dead_particles

    dt = 5.3e-12 #step approximately equates to 1cm per frame

    for p in active_particles[:]:
        p.update_active(dt,dims) 
        if p.absorbed:
            active_particles.remove(p)
            fading_particles.append(p)

    for p in fading_particles[:]:
        p.update_fading()
        if p.dead_ring_size == 0:
            fading_particles.remove(p)
            dead_particles.append(p)

    sim_time = frame * dt
    timer_text.set_text(f"Time: {sim_time*1e9:.2f} ns")
    alive_text.set_text(f"Alive: {len(active_particles)}")
    pmt1_text.set_text(f"CH1: {pmt_hits[0]}")
    pmt2_text.set_text(f"CH2: {pmt_hits[1]}")

    artists = []
    for p in active_particles + fading_particles:
        artists.extend([p.dot, p.trail_line])

    return artists

def plotter(dims,pmts):
    fig,ax = plt.subplots()
    ax.set_xlim(-10, dims[0]+10)
    ax.set_ylim(-10, dims[1]+20)
    ax.set_aspect('equal')

    box = plt.Rectangle((0,0),dims[0],dims[1],fill=False)
    ax.add_patch(box)

    for p in pmts:
        width = p['xmax'] - p['xmin']
        height = p['ymax'] - p['ymin']
        ax.add_patch(plt.Rectangle((p['xmin'], p['ymin']),width,height,color='green'))
        ax.text((p['xmin']+p['xmax'])/2,(p['ymin']+p['ymax'])/2,str(p['channel']),fontsize=12,fontweight='bold',ha='center',va='center')

    timer_text = ax.text(0.02, 0.98, "", transform=ax.transAxes, fontsize=10, va="top", ha="left")
    alive_text = ax.text(0.02, 0.92, "", transform=ax.transAxes, fontsize=10, va="top", ha="left")
    pmt1_text = ax.text(0.4, 0.98, "", transform=ax.transAxes,fontsize=10, va="top", ha="left", color="red") 
    pmt2_text = ax.text(0.54, 0.98, "", transform=ax.transAxes,fontsize=10, va="top", ha="left", color="red") 

    return fig,ax,timer_text,alive_text,pmt1_text,pmt2_text

#randomized functions

def particlegen(n,pos,ax):
    particles = []
    for _ in range(n):
        theta = np.random.uniform(0,2*np.pi)
        velocity = v*np.array([np.cos(theta),np.sin(theta)]) 
        particles.append(Particle(pos,velocity,ax))
    return particles


"""SETUP"""

#PMTs
dims = np.array([110,76]) 
pmts = [
    { 'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': 76,   'ymax': 81,  'channel': 1 },
    { 'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': -5,   'ymax': 0,   'channel': 2 },
]

#plotting things
trail_length=1 

#physics things
TIR = np.radians(39)
attenuation = 380


"""OUTPUT"""


#output folders
output_folder = "anglescans"
gif_folder = os.path.join(output_folder, "gifs")
os.makedirs(gif_folder, exist_ok=True)

#CSV file
csv_file = os.path.join(output_folder, "anglescans.csv")
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["run", "x", "y", "pmt1_hits", "pmt2_hits"])


"""Scanning Loop"""


ys = [38.5] #y values to scan
xs = np.arange(2.5, 7.5, 5)
run = 1

for y in ys:
    for x in xs:
        start_time = time.time()

        #reset globals
        pmt_hits = [0,0,0,0]
        fading_particles = []
        dead_particles = []


        print(f"Starting Run{run} ({x},{y})")
        fig, ax, timer_text, alive_text, pmt1_text, pmt2_text = plotter(dims, pmts)

        active_particles = particlegen(5000,[x,y],ax)

        animation = FuncAnimation(fig,animate,interval=10,frames=500,blit=True,init_func=init)

        gif_path = os.path.join(gif_folder, f"run{run}_x{x}_y{y}.gif")
        animation.save(gif_path, writer='pillow', fps=10)

        #log results
        with open(csv_file, "a", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([run, x, y, pmt_hits[0], pmt_hits[1]])

        end_time = time.time()
        print(f"Run {run} (x={x}, y={y}) completed in {end_time-start_time:.2f} s")

        run += 1