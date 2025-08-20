This is a bit of documentation to help you navigate the python I used for my work in the milliQan experiment over the summer of 2025.
It will include: 
	Some clarity on how each portion of the code operates 
	Simple instructions for some alterations that may be necessary
	Common bugs/issues you may run into and how to fix them.

Notebook for Data Analysis:

	Note: Not really any bugs or instructions, more just running through what I did. most of this works and updates itself with new data
	
	Importing data:

		The first cell in this notebook is used to import data from a txt file in the cms computers and converting it to a pandas dataframe
		When manually sshing into cms1, you use ssh {username}@cms1.physics.ucsb.edu, then input your password. this is what the username and password is referring to
		After sshing in, we use cat to grab the file, send it to the 'output' variable, strip it of lines that may cause any errors, and load it into a pandas data frame
		This data frame is saved as df_base
		The data frame is organized into 9 columns: [Run#, x, y, HVSlab, HVBar, MPV3, dMPV3, MPV4, and dMPV4]
		MPV are the recorded peak position of a pulse area histogram, and dMPV is the error in that recording

	First process:
		
		The data frame that we have at this point has a number of gaps where, for one reason or another, there is no actual data to analyze
		We eliminate those gaps in the second cell, where we grab the last 4 columns and check if they have any results in them. If not, we remove them

	Ch 3 analysis:
		
		Power Law Fit:

			Before we can do any real analysis, we must ensure that our data does not have any dependence on HV
			To do this we fit it to a relationship between MPV3 and HVSlab, then extrapolate that fit out to 1450V, effectively removing HV dependence in the runs
			First I grab the data corresponding to a set of runs that have the same position, but different HVSlabs
			This is done by grouping lines of data by their position (x,y), and then asking if a group has more than one unique HVSlab value
			Every line that satisfies this into the data frame labeled 'multi'
			The fit we are applying looks like ln(MPV3) = ln(a) + b*ln(HVSlab)
			It's easiest to perform this fit by logging all the data we put in, so we do that
		FIGURE THIS OUT FILL IT IN

		Note: I'm not sure how error should be processed through this, so I skipped over it

		2-d Histograms and Interpolation:

			From this point, now that we have our fitted data, we can use it basically however we want
			First we extract from the data frame the x values, y values, and extrapolated MPV3
			We then use these values to create a 2D histogram with matplotlib
			We do not have the limits of the plot defined, so if you use a smaller data set it will not fill out to the full slab size
			Arguments of note:
				c is the z axis in this 2D hist, its what the coloring represents
				s is the size of the dots on the histogram
				gca().invert_yaxis() inverts the y axis. note that this must come before applying labels and a title or else it will not work

		
			Now I want to connect the dots
			To do this, we will interpolate
			First we log our data to shrink its effective scale. This helps scipy make accurate interpolations
			Then we must set up a mesh, defining little boxes that we will assign interpolated values
			we now use scipy's griddata to estimate the values between are known points.
			cubic is the most accurate of the 3 interpolation methods
			We then convert the interpolated data back into a normal scale
			This Zgrid is then plotted onto a color mesh along with the actual data, using the same the same methods.

		Angle Dependence plots:

			For these plots, I want data from scans across the slab at single y values
			To achieve this, we use masks
			We assign a mask to a criteria (in this case y = 65)
			Then we draw out the data we want that 'wears' the mask
			Now when calculating the angle of the position relative to the PMT, there are a couple quirks that arise due to the method of how I recorded position:

				I recorded positions by taking the inside measurement from the edge of the slab to the bar PMTs
				This means I recorded the position of a corner of the selection area, as opposed to its center
				To account for this in the calculation, I add 2.5 to both x and y so that I can use the center of the face of the PMT as a reference
				i.e. for the point (20,65), i use the point (22.5,67.5) in my angle calculation.
				The center point of the faces of the PMTs are located at (30,76) for CH3 and (30,0) for CH4

			Now we just plot MPV3 against angle in degrees
			If you want to add more scans, just copy paste the masking method, angle calculation, and plotting, replacing any instance of the old y value with what you now want
			If scans you added do not show up/throws you some error that it doesn't exist, you probably haven't updated the original df yet. Run the notebook top to bottom

	CH4 Analysis:

		Much of the CH4 Analysis is identical to that of CH3, however one main difference is a reduced data set
		Early runs in the data set do not have accurate MPV4 values, so we remove them from the data frame
		If for one reason or another you wish to use all the data in the df, simply comment out the line: df_ch4 = df[df["Run"] >= 65903] from the first cell of the section

		Power Law Fitting:

			We use an identical fitting method, just swapping 3s for 4s in the code
		Histograms and Interpolation:

			Similarly an identical histogram and interpolation method

		Angle dependence plots:
	
			I use the same addition to my x and y to calculate the angles for those plots, however now I use CH4 located at (30,0)

	Simulation Analysis:
		Input data:
	
			Using the automodel.py script writes a .csv that we read in the first cell
			It is formatted [run#,x,y,pmt1_hits,pmt2_hits]
			Note that run# here is not related to the run# from real data. Here it is simply for bug fixing purposes

		histograms and Interpolation:
			We extract the data we want
			We use a very similar 2D histogram plotting method as before
		
		Angle Dependence Plots:
			Here we define the PMT positions so that if we run a simulation with a different slab geometry, adjusting the code is simple
			The angler function effectively does what we were doing manually before
			To avoid clutter for if you want to plot a lot of scans together, this is organized into a for loop
			Define the number of columns you want and the for loop will grab the scans at each y value and plot them. together
			When we combine the real data with the simulation data, you first have to define which scans you want to grab from the simulation. you do not have to have the exact same y values or even the same number of y values that exist in the actual data. It will only plot what it has. The 3 for loops plot the sim data, then the real data, then sets up some plot settings.


Common Bugs:
	Not all the modules are imported in every cell, nor is all the data. It works if you go top to bottom, so I would suggest when first opening this notebook to just run every single cell and then begin to make edits. 
	The other common bug you will run into is when first importing the data from cms2, it might say that it failed because there weren't and rows or something. Sometimes it takes longer than we expect for the ssh-ing to go through and the cell accidentally looks for a file that doesn't exist and will provide nothing. Just run this cell again, and eventually it will work. 
	

			
		