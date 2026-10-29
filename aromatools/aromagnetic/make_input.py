#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

GEN_INPUT = """
# ================================
#        AroMagnetic INPUT
# ================================

# --------------------------------
#      General Information
# --------------------------------

NAME_MOLECULE             name_molecule
PROGRAM                   ORCA                      # ORCA or GAUSSIAN
KEYWORD_LINE              B3LYP def2-TZVP NMR 
PROCS                     8
MEMORY_GB                 16
CHARGE                    0
MULTIPLICITY              1
MOLECULAR_COORDINATES     mol.xyz
PSEUDO_PI                 NO                        # YES or NO

# --------------------------------
#      GRID (Grid generation)
# --------------------------------

BEGIN_GRID
 
GA_MAX            200                               # Maximum number of ghost atoms per input
EXTERNAL_FIELD    0.0 0.0 1.0                       # Bx By Bz
NUM_POINTS        50 50 50                          # Nx Ny Nz
GRID_MARGIN       2.0                               # Å
VTK_FILES         YES                               # YES or NO

END_GRID

# --------------------------------
#        VERTICAL_PROFILE
# --------------------------------

BEGIN_VP

SCALE       Linear                                  # Linear or Logarithmic
START       0.0
DELTA_Z     0.1
NUM_GA      100
# END       5.0                                     # Only for Logarithmic scale

# ================= OPTIONS =================
RCS         YES                                     # YES or NO
INICS       YES                                     # YES or NO
NICS_SCAN   YES                                     # YES or NO (with z-profile)

END_VP

# --------------------------------
#       HORIZONTAL_PROFILE
# --------------------------------

BEGIN_HP

POINTS           81
MARGIN           2.0
HEIGHT           0.0 1.0 1.7
DELTA_X          0.25
MAKE_AXIS_Y      YES                                # YES or NO Just for NICS_SCAN_XY
DIRECTION        FORWARD                            # FORWARD or BACKWARD Just for PREDI_XY

# ================= OPTIONS =================
NICS_SCAN_XY     YES                                # YES or NO
PREDI_XY         YES                                # YES or NO

END_HP

# --------------------------------
#        Orbital Partition
# --------------------------------

BEGIN_ORBITALS

MO_core     1 2 3
MO_sigma    4 5 6
MO_pi       7 8 9
MO_others

END_ORBITALS


# --------------------------------
# ORCA Block (Optional)
# --------------------------------

BEGIN_ORCA_BLOCK

%pal
  nprocs 8
end

%maxcore 2000

END_ORCA_BLOCK

"""

def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else "INPUT.txt"
    with open(filename, "w", encoding="utf-8") as f:
        f.write(GEN_INPUT)
    print("Done!!")
    print("Now, you should complete with your specifications the INPUT.txt file!")
  
if __name__ == "__main__":
    main()

