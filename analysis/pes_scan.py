import numpy as np
import matplotlib.pyplot as plt

# plot the energies for a bunch of relaxed defect configurations (generated via the interstitial.py code e.g.)
# with respect to some configuration coordinate Q (e.g. RMS distance of halides to central Pb ion, or angle of rotation around Pb)
# (avg disp between pristine/relaxed geometries doesn't seem to give much meaningful correlation,
# need coordinates which are sig diff between configs)
# this will provide a 1D scan of the PES of the defect

# could also do a 2D scan by giving 2 configuration coordinates (e.g. distance & angle? or coordinates along 2 lattice vectors),
# would need a 2D colourmap or a 3D surface plot
# OR could also give all 3 cartesian/direct coords of defect ion wrt pristine geometry,
# then CC plot would be a 3D surface plot w colourmap 

# need to read the total energy and structure outputs from VASP
# so need to parse the OUTCAR and CONTCAR files...might actually be easier to grep with some shell code?

# if i want to use RMS dist of all halides surrounding a Pb ion, need to identify the defect ion, then nearest Pb ion,
# then identify all 6/7 closest halides
# already have code for sorting halides by distance around given Pb atom, steal from interstitial.py

