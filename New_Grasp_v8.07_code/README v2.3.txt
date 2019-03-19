Bayesian program

Installation instructions:

Extract .zip file to a directory.

Before running it the first time, go to ./minFunc_2012 and run mexAll.m
This will compile some stuff to make the minFunc routine faster.

Instructions for use:

Use Grasp v6.79 or later.
Start Grasp from the installation directory using the commands:
'grasp_startup_6.79;grasp'

You may need to edit grasp_startup_6.70.m if you have multiple versions of Grasp or are using a non-standard directory.

Load in the data as usual (unfortunately you will have to start a new project, it won't work with existing projects).

Find a value for the rocking width in the usual way for a particular spot with coordinates [x0 y0]

Create a section in run_OO.m below line 44 in the form:
	
	case 'Ca doped YBCO 16.4T phi'
            input_index = 2; %data location (number)
            output_index = 2; %where to put result
            eta0 = 0.79139;  %rocking width for spot in degrees
            spot = [181 105];  %Spot coordinates
            sanoffset =0;  %san misalignment
            phioffset =0;  %phi misalignment


Change "inputs = {'BFAP 12T'}" on line 22 to point to this case.  You can use more than one for combined results.

Things to try:

fit = 1  % tries to fit rocking width and/or misalignment only for data in:
sectors = 0, 'sector_mask' or'sector_boxes' (line 34), meaning all pixels, sector, or sector boxes respectively.  Quite important for fitting, it will give a bad answer if you give it too much irrelevant stuff.

Fitmethod - try 'check' (fminunc if optimization toolbox installed) or 'minFunc'. The latter may work even if fminunc fails, and has lots of options.

parameters e.g. offsets can be fixed, for example after using Nb to align.

"inputs = {'data1', 'data2'}" on line 22 will use posterior of 'data1' as prior of 'data2'.  This can be used to combine san and phi rocks.  In that case eta0 in 'data2' will be ignored, and the rocking width (in q) given (or fitted) from data1 will be used.

If fitting is used, the results will be output to "results.txt" in the format:
'inputs' rock_width(in q) err rock_width(degrees) err sanoffset err phioffset err


A box sum for a given spot should give an integrated intensity in q, i.e. the angular integrated intensity with the lorentz factor already included and multiplied by the value of mod(q).  I have not yet verified whether this is correct!  I'd like to know very much if corresponds to the number from more traditional methods!

