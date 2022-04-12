from trojan import *

'''
The parameters chosen in levels are the min, max and spacings of contour lines.
If this is not included, the contours will sink into the centre (where the Sun is) and look like a solid block

'''

planet = Two_Body_System(0.01, 1)

planet.planar_plot(-1.5,1.5,-1.5,1.5)
