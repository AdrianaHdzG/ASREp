import sys
from pathlib import Path
import numpy as np

import os

ASREp_dir = "/Users/jinyanzhao/Desktop/tools/ASREp"
if not os.path.exists(ASREp_dir):
    sys.exit('Error: Please set the ASREp_dir variable to the absolute path of the ASREp directory')
sys.path.append(str(ASREp_dir))

import platform
print(f"Platform processor: {platform.processor()}")

import ASREp
import ASREp.ground_deformation_models as gdm
from params import *

# Define beam and tunnel parameters 
# (vlt and Kt, Eb, nu_b, E_s are uncertain params defined in the params.py)
beamX = np.linspace(-10, 10, 41)
z_t = 8.0  # Tunnel depth
D_t = 6.0   # Tunnel diameter

EoverG = 2+2*nu_b

dfoot = 9
bfoot = 0.5
qfoot = 3.2 * 10 * 1000
nus = 0.5
mu_int = np.tan(30*np.pi/180)

(vertical_displacement_ground_building, horizontal_displacement_ground_building
     ) = gdm.ground_disp_Mair_1993(
        beamX,
        z_t,
        Kt,
        D_t,
        vlt
    )


beamY = np.zeros_like(beamX)
beamZ = np.zeros_like(beamX)

solver = 'EL'  # or 'EP', depending on the desired solver

if solver == 'EL':
    model = ASREp.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                     beamZ, dfoot, bfoot,
                                     solver='elastic')
    model.set_beam_properties(Eb, EoverG, qfoot)
    model.set_soil_properties(Es, nus, mu_int)
elif solver == 'EP':
    model = ASREp.ASRE_Timoshenko_model(beamX.size, beamX, beamY,
                                     beamZ, dfoot, bfoot)
    model.set_beam_properties(Eb, EoverG, qfoot)
    model.set_soil_properties(Es, nus, mu_int)
else:
    sys.exit('Error: string solver must be EL or EP')

model.run_model(
    horizontal_displacement_ground_building, np.zeros_like(horizontal_displacement_ground_building),
    vertical_displacement_ground_building, 'strain+disp'
)

max_tensile_strain = model.eps_t.max()

with open('results.out', 'w') as f:
    f.write('{:.60g}'.format(max_tensile_strain))