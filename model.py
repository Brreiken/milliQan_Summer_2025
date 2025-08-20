import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from collections import deque
import time

matplotlib.use("Agg")

"""Globals"""

#constants
c = 299792458 #speed of light
n = 1.58 #refractive index of ej200 plastic scintillator
v = c/n

#empty sets
pmt_hits=[0,0,0,0]
fading_particles = []
dead_particles = []

"""Functions"""

class Particle: #particle identity stuff
    def __init__(self,position,velocity, ax):
        self.position = np.array(position,dtype=float) #defines the particles position
        self.velocity = np.array(velocity,dtype=float) #now its velocity
        self.dot, = ax.plot([],[],'o', markersize=4,color='blue') #and the visualization of the particle

        self.trail_length = trail_length #defining the trail length using a global defined later
        self.positions = deque(maxlen=trail_length) #this makes the positions array limit itself to only ever be so long. i.e. for a maxlen=50, the array will fill positions[0]...positions[49], then positions[1]...positions[50], and so on
        self.trail_line, = ax.plot([],[],'-',lw=1,color='blue') #setting up the trail line
        self.distance_traveled = 0.0  # in cm

        self.dead_ring_size = 4  # starting marker size when absorbed
        self.dead_ring_fade_rate = 0.125  # how much it shrinks per frame

        self.absorbed = False 
        self.counted = False

    def pmt_check(self):
        for p in pmts:
            if (p['xmin'] <= self.position[0] <= p['xmax'] and p['ymin'] <= self.position[1] <= p['ymax']): #if within the dimentions of a PMT
                self.counted = True
                pmt_hits[p['channel'] - 1] += 1 #add the count to its respective PMT
                #print(f"Channel {p['channel']} hits: {pmt_hits[p['channel'] - 1]}") #Still for now just shows up in the terminal, i suggest commenting out when doing summing

                #sumtop = pmt_hits[0]+pmt_hits[1] 
                #sumbot = pmt_hits[2]+pmt_hits[3]
                #print(f"Top sum hits: {sumtop}")
                #print(f"Bot sum hits: {sumbot}") #These 4 lines are useful for summing with 4 PMTs. If you want it, just uncomment the lines

                break

    def update_active(self,dt,bounds): #update every frame
        if self.absorbed:
            return

        dx = self.velocity*dt*100 #*100 is to convert velocity from m/s to cm/s
        self.position += dx

        survival_prob = np.exp(-np.linalg.norm(dx)/attenuation) 
        if np.random.rand() > survival_prob: #rolls the dice
            self.absorbed=True
            return
        
        for i in range(2): #2 dimentions for the walls

            #right wall = die, just comment out to turn it off
            if i == 0 and self.position[0] >= bounds[0]:
                self.absorbed = True
                return

            if not 0 < self.position[i] < bounds[i]: #If not in box
                self.velocity[i] *= 1 + np.random.uniform(-0.04,0.04)
                self.velocity[i-1] *= 1 + np.random.uniform(-0.04,0.04) #these perturbations are an effective equivalent to defining surface roughness
                angle = np.arctan2(abs(self.velocity[1-i]),abs(self.velocity[i])) #defines the angle between the particle and the normal to the edges of the box
                if angle < TIR: #check total internal reflection
                    self.pmt_check() #if absorbed, is it in the pmt?
                    self.absorbed = True
                    return
                else:
                    self.position[i] = np.clip(self.position[i],0,bounds[i]) #make sure that the particle gets put back in the box if the jump from frame to frame would send it out.
                    self.velocity[i] *= -1 #bounce
                    self.velocity *= v/np.linalg.norm(self.velocity) #make sure the pertubation doesn't get us faster than v

        self.positions.append(self.position.copy()) #updating positions for the tail
        x,y = zip(*self.positions)
        self.trail_line.set_data(x,y)
        self.dot.set_data([self.position[0]], [self.position[1]])


    def update_fading(self):
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
    
#animation based functions

def init(): #sets up the initial frame of the animation
    artists = [] #dots on the animation are 'artists' 
    for p in active_particles + fading_particles:
        p.dot.set_data([p.position[0]], [p.position[1]])
        p.trail_line.set_data([], [])
        artists.extend([p.dot, p.trail_line])
    return artists

def animate(frame): #frame here serves to effectively iterate this function. it will run each frame
    global lastUpdate, active_particles, fading_particles, dead_particles
    #dt = time.time() - lastUpdate #dt is the real time between 2 frames
    #lastUpdate=time.time()

    
    dt = 5.3e-12 #This works out so that each frame is approximately a 1 cm step

    for p in active_particles[:]:
        p.update_active(dt,dims) #run the update each substep
        if p.absorbed:
            active_particles.remove(p)
            fading_particles.append(p)

    for p in fading_particles[:]:
        p.update_fading()
        if p.dead_ring_size == 0:
            fading_particles.remove(p)
            dead_particles.append(p)

    sim_time = frame * dt  # in seconds
    timer_text.set_text(f"Time: {sim_time*1e9:.2f} ns")  # ns
    alive_text.set_text(f"Alive: {len(active_particles)}")
    pmt1_text.set_text(f"CH1: {pmt_hits[0]}")
    pmt2_text.set_text(f"CH2: {pmt_hits[1]}")
    #pmt3_text.set_text(f"CH3: {pmt_hits[2]}") 
    #pmt4_text.set_text(f"CH4: {pmt_hits[3]}") #uncomment just these 2 lines if you only want 4 pmts
    #sumtop_text.set_text(f"Top sum: {pmt_hits[0]+pmt_hits[1]}")
    #sumbot_text.set_text(f"Bot sum: {pmt_hits[2]+pmt_hits[3]}") #uncomment all 4 lines if you want 4 pmts and the summing effect

    artists = []
    for p in active_particles + fading_particles:
        artists.extend([p.dot, p.trail_line])

    return artists

def plotter(dims,pmts): #set up the plot for the animation
    fig,ax = plt.subplots()
    ax.set_xlim(-10, dims[0]+10)
    ax.set_ylim(-10, dims[1]+20)
    ax.set_aspect('equal')

    box = plt.Rectangle((0,0),dims[0],dims[1],fill=False) #our slab in 2-d :)
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
    #pmt3_text = ax.text(0.68, 0.98, "", transform=ax.transAxes,fontsize=10, va="top", ha="left", color="blue")  
    #pmt4_text = ax.text(0.82, 0.98, "", transform=ax.transAxes,fontsize=10, va="top", ha="left", color="blue")  
    #sumtop_text = ax.text(0.4, 0.92, "", transform=ax.transAxes, fontsize=10, va="top", ha="left", color="red")
    #sumbot_text = ax.text(0.68, 0.92, "", transform=ax.transAxes, fontsize=10, va="top", ha="left", color="blue")

    return fig,ax,timer_text,alive_text,pmt1_text,pmt2_text,#pmt3_text,pmt4_text,sumtop_text,sumbot_text

#randomized functions

def particlegen(n,pos,ax):
    particles = []
    for _ in range(n):
        theta = np.random.uniform(0,2*np.pi) #generates random angles
        velocity = v*np.array([np.cos(theta),np.sin(theta)]) #the particles are all at speed v, each having one of these random angles
        particles.append(Particle(pos,velocity,ax)) #add them to the list
    return particles

"""main"""

#Slab stuff
dims = np.array([110,75]) #lxw, (x,y)
pmts = [
    { 'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': 75,   'ymax': 80,  'channel': 1 },
    #{ 'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': -5,   'ymax': 0,   'channel': 3 },
    #{ 'xmin': 60-2.5, 'xmax': 60+2.5, 'ymin': 75,   'ymax': 80,  'channel': 2 }, 
    #{ 'xmin': 60-2.5, 'xmax': 60+2.5, 'ymin': -5,   'ymax': 0,   'channel': 4 },
    { 'xmin': 30-2.5, 'xmax': 30+2.5, 'ymin': -5,   'ymax': 0,   'channel': 2 }, #if you want 4 pmts instead of 2, uncomment the 3 lines above this and comment this line
]

#plotting stuff
trail_length=1 #shorter tail might be necessary for efficiency

#physics variables
TIR = np.radians(39) #approximate calculated angle of total internal reflection between scintillator and air
attenuation = 380 #EJ-200 attenuation length in cm for peak wavelength ~425nm

#actual doing stuff
#fig, ax, timer_text, alive_text, pmt1_text, pmt2_text, pmt3_text, pmt4_text, sumtop_text, sumbot_text = plotter(dims, pmts)#make the plot
fig, ax, timer_text, alive_text, pmt1_text, pmt2_text, = plotter(dims, pmts)#2 pmt plot

active_particles = particlegen(500,[35,40],ax) #generate the particles

start_time = time.time()
lastUpdate=time.time()
animation = FuncAnimation(fig,animate,interval=10,frames=500,blit=True,init_func=init)
animation.save("particle_sim3.gif", writer='pillow', fps=10)
#plt.show()
end_time = time.time()
print(f"{end_time-start_time:.2f}")
