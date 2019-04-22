# femtography
Source codes for Femtography project
To run it use:
./bin/PARTONS_example MODE OUTPUT_FILE

The available MODEs are:

GPDGK16Numerical
Evaluate GK16 model for H_u in points of (x, xi, -t) for a fixed scale

GPDMMS13
As before but for MMS13 model

CFF
Evaluate DVCS CFF H for LO and GK16 model in points of (xi, -t, Q2)

OBS_ALU
Evaluate DVCS ALU asymmetry  for LO CFFs and GK16 model in points of (xB, -t, phi) for a fixed scale and beam energy

OBS_CS
Evaluate DVCS unpolarized cross-section  for LO CFFs and GK16 model in points of (xB, phi, subProc) for a fixed t, scale and beam energy, where
subProc = 0 BH
subProc = 1 INT
subProc = 2 DVCS

You will easily spot the place in the code where one can change ranges, fixed values and number of steps. 
For OBS_ALU and OBS_CS it is checked if a given kinematic point is valid, e.g. it is checked if 0 < y < 1
3+1 D distributions are produced, which should work for the visualization.
If you want to speed up computations, you can modify computation.nb.processor in bin/partons.properties. The multi-threading will be taken care of automatically.
