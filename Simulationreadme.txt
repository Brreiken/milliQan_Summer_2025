This is a bit of documentation to help you navigate the python I used for my work in the milliQan experiment over the summer of 2025.
It will include: 
	Some clarity on how portions of the code operate
	Simple instructions for some alterations that may be necessary
	Common bugs/issues you may run into and how to fix them.

Note: I will not include specifics on any plotting done in the various scripts and notebook as for the most part its personal preference of how you use matplotlib and its tools.

model.py
	-This is the first, simple model that I made, and all others at a basic level follow its processes.
	-At a fundamental level, it operates via the animate function, which runs for a number of frames and each frame sends a series of checks and updates
	
	-Particle Class
		-init function
			does not explicitly get called, but instead automatically runs when defining the particle class. Sets up some basic phyics and plotting properties for the photons
		-pmt_check
			Once our particle hits a wall, we ask whether its inside the dimensions of any of the PMTs. 
			This function currently contains some commented out lines that  describe the summing effect we have on the in situ slab detector. There are similar commented out lines in other functions that are required for these to work.
		-update_active
			This first makes sure the photon is still alive, then steps the photon forward, defines the distance traveled in that step, and calculates the chance for the photon to survive attenuation. We say P = e^(-x/A), where A is the attenuation length. We then quickly roll the dice on whether or not the photon will be attenuated. After this we have the wall checks. First I define whether or not one end of the slab is going to absorb all of the photons, like we assume for the Broida detector, then I just ask if after taking our step forward from earlier if we're still inside the box. If we are not, we check for total internal reflection, where we use a roundabout method for surface roughness. Here we perturbate the velocity, and thus the angle, of the traveling photon right before we check for total internal reflection. This is an effective equivalent to making the surface of the scintillator rough, and simplifies the photon velocity change we would see from reflecting off of an angled surface. If the photon is within the angle of total internal reflection, it will be absorbed into the wall, and run our PMT check. If the photon is outside the angle of total internal reflection, it will be pushed back into the slab, reflect using its perturbed velocity, and renormalized to ensure that our photons are not exceeding the speed of light. We finally then update the position of the photon.
			The single intended alteration that exists in this portion of the code is you can comment out the lines saying one wall will absorb all photons. It is noted in the code itself.
			Attenuation and surface roughness are isolated, if you remove them, the code will operate fine without them. Total internal reflection, however is a property that you cant just comment out to remove. Instead if you wish to remove it, change the TIR variable defined later to 0 if you wish for the photons to always reflect or anything greater than 90 if you wish for the photons to always absorb.
			All the physics is contained in this function, so if the code runs, but the photons do not behave how you want them to, you likely need to make changes here.
		-update_fading
			This is a function that defines a visual effect for the photons in the plot. It makes it so that absorbed photons will fade and not clutter the plots.
	
	-Animation functions
		-init
			This initialization function is for graphing purposes, it sets up the frame and the photons initial positions before any of the loops
		-animate
			This function serves as a core loop for the entire script. We define the time step between each frame, then run the updates on all active and fading photons.
			Adjusting dt will change how far the photons move each frame. This can reduce run time, but can also reduce physical accuracy
			This function also contains 4 lines that add counters for 4 PMTs and for summing into the graphic. Uncomment them if you want 4 PMTs and summing

		-plotter
			This function actually makes the plots that we use.
			This function also contains 4 lines that add counters for 4 PMTs and for summing into the graphic. uncomment them if you want 4 PMTs and summing

	-Randomizing Functions
			-particlegen
				This is currently the only randomizing function we have. It simply creates n photons at a specified initial position and provides them with velocities in random directions

	-MAIN
		This main section contains all of the defined variables not in the particle class and the start of all the loops.
		The pmts variable contains 3 lines that define a setup for 4 PMTs and one line that only works for a 2 PMT setup. To swap to a 4 PMT setup, uncomment the 3 lines and comment the one.
		A number of variables here can be freely changed to alter physical effects or positions of objects. The majority are straightforward. The only one that is not really is active_particles. Here you can change the number of photons produced and where by changing the first two arguments input to particlegen. The position must be input in an array [x,y].



automodel.py
	This script is nearly identical to model.py, however instead of simulating for a single position, it ends in a for loop that allows you to simulate over many positions
	This script does not include many of the comments from model.py. In particular if you wish to change to a 4 PMT geometry, refer to model.py
	-output
		we define an output folder for everything from your scan to go into. It will go into the same directory as where automodel.py is located. Inside that folder we make another one to put all the gifs that are made. We then make a .csv that will contain all of the actual counting data that we want to record.
		Altering names is recommended between scans. If you alter the .csv here you must also change the data recording in the for loop to match it
	-final for loop
		First we define the y values you wish to scan across, then the x values you wish to scan at. Then, in the for loop, we run the particlegen and plotter functions as in model.py, but iterated over each y and then x you specified. After each run, it opens the .csv and adds a new line with all of the data that you wish to record.
		The positions that you specify must be provided in arrays for this to work. If you want to scan across the slab holding x constant, simply flip the for y in ys: and for x in xs: similarly to model.py if you wish to alter any of the variables at the end, they're all pretty self explanatory with the number of photons generated still being contained in the particlegen input.
		Run# is something I used strictly for bug testing. If you wish to remove it, simply remove 'run' from writerow([...]) when making the .csv and in the for loop

3dmodel.py
	This script is very similar to model.py. Instead of all the 2D physics and plotting that exists there, though, there is now 3D physics. As I will not be detailing the specifics of the plotting, I will only go through the alterations in the physics between this script and model.py.
	
	-pmt_check
		This function now checks all 3 dimensions to see if the photon is inside the PMT. Otherwise operates identically to model.py
	-update_active
		We now have 3 dimensions of walls. The check to see if the photon is inside now just includes the third dimension. The main difference between model.py and 3dmodel.py exists here. Surface roughness and total internal reflection are now handled differently. Now what we do is define a normal vector with an inward facing normal based on the wall that we hit. We then perturb this normal vector to simulate the varying surface angles. When then define the components of the photon's velocity parallel to the inward facing normal and parallel to the wall, and calculate are angle with theta = arctan(v_parallel/v_normal). We then check this angle against the total internal reflction angle and do the absorption/reflection responses. If the photon reflects, it must now reflect off of the angled surface, which is done by using something like v' = v - 2*normal*perturbed_normal.
		This function now contains a new variable in it called pert_sigma, which is the standard deviation of the distribution of perturbations we get in this function
	-particlegen
		this function operates in a similar manner, however now we use spherical coordinates to generate random directions we randomize an azimuthal angle, then randomize the cosine of the polar angle and calculate the corresponding sine. We use these 3 to transform a direction from spherical to cartesian, and then apply our velocity
	-main
		now all of our dimension based variables have 3 dimensions to worry about. They all vary the same way. pert_sigma is located here.

	Generally speaking, all plotting tools were adjusted so that we create 3 2D frames, one for each plane.

	

Common shared bugs:
	Between all three of these scripts, the single most common issue you can run into is changing only some of the necessary components for a new geometry but not all. 
	For example if you comment/uncomment some of the lines to change from a 2 PMT geometry to a 4 PMT geometry you will likely get thrown an error somewhere of a variable not being defined, whether that's because you commented out its definition or because you forgot to return it. Follow the error message to the line it asks for and it's typically a really simple fix. 
	Other geometry issues may be changing the shape of the slab, but not changing the positions of the PMTs. You may have a run that gives you 0 PMT counts that you can't explain away with physics. Check the gif you made and see where the PMTs are relative to the slab. if they do not at the very least share a wall, the photons cannot get into the PMTs and cannot get counted. Adjust the 'pmts' dictionary accordingly
	Other than these two, 99% of bugs you will run into will come as a result of you adding a physical effect and it having a poor interaction with some other part of the code.

	Note: In all these simulation scrips, the physics is all tied into the animation function, which also wants to make a plot. If you wish to avoid making the plot and save on some processing time, you must manually define a new for loop that will run at the very least run update_active in a similar form.

	