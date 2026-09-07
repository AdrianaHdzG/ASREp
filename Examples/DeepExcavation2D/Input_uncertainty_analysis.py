# This script defines the uncertainty variables used in quoFEM

# MATERIAL PARAMETERS

EA      = 25.9                        # [GN]
GAs     = 5.95                        # [GN]
EI      = 2777                        # [GN*m^2]

# MATERIAL PARAMETERS
poissons_ratio      = 0.3           # [-] For the building materials

# Wall parameters
avg_wall_disp                 = 0.0497     # [%] Beta
retaining_wall_depth          = 14      # [m] Depth of the retaining wall from surface to rigid layer

# SOIL INITIAL PARAMETERS

Es                          = 175        # [MN/m2] (COV = 14-68% according to Zhao and DeJong, 2023. The lower limit is used --> standard deviation = 0.14) 175 MN/m2 is for ML from Structural Design and 3MpA for the Fill. Foundation at +2.7 and top of ML at +2. 
