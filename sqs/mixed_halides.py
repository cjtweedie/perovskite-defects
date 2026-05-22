import numpy as np
import matplotlib.pyplot as plt
from math import factorial 
from ase import Atoms
from ase.spacegroup import Spacegroup
from ase.visualize import view
from ase.build import make_supercell
from ase.spacegroup import crystal
from ase.io import read
from ase.io.vasp import write_vasp
from icet import ClusterSpace
from icet.tools.structure_generation import (generate_sqs, 
                                             generate_sqs_from_supercells, 
                                             generate_sqs_by_enumeration, 
                                             generate_target_structure)
from icet.input_output.logging_tools import set_log_config
set_log_config(level='INFO')
from generate_sqs import generate_sqs

# x = 1/3 -> CsPbI2Br
# alloy lattice consts interpolated via Vegard's law from pristine CsPbI, CsPbBr DFT-relaxed cells 
pure_cubic_13 = crystal(('Cs', 'Pb', 'I'), basis=[(0.5,0.5,0.5), (0.0,0.0,0.0), (0.0,0.0,0.5)], spacegroup=221, cellpar=[6.16,6.16,6.16,90,90,90])
pure_tet_13 = crystal(('Cs', 'Pb', 'I', 'I'), basis=[(0.5,0.0,0.5), (0.0,0.0,0.0), (0.0,0.0,0.5), (-0.212,0.288,0.0)], spacegroup=127, cellpar=[8.49,8.49,6.27,90,90,90])
# fix up ortho, basis is incorrect I think? -> doesn't interpret the space group properly, since in non-standard setting? but wyckoff positions also wrong idk
#pure_ortho_13 = crystal(('Cs', 'Pb', 'I', 'I'), basis=[(0.25,0.9835,0.4451), (0.5,0.0,0.0), (0.5306,0.1973,0.2993), (0.75,0.9369,0.9983)], spacegroup=62, cellpar=[8.32,8.77,12.24,90,90,90])
pure_ortho_13 = read("ortho_111.cif")
#print(pure_ortho_13)
#view(pure_ortho_13)

# order of elements for CsPb(I_xBr_1-x)3
CsPbIBr_order = {
    'Cs': 0,
    'Pb': 1,
    'I': 2,
    'Br': 3,
    }

# x = 2/3 -> CsPbIBr2
# alloy lattice consts interpolated via Vegard's law from pristine CsPbI, CsPbBr DFT-relaxed cells
#pure_cubic_23 = crystal(('Cs', 'Pb', 'I'), basis=[(0.5,0.5,0.5), (0.0,0.0,0.0), (0.0,0.0,0.5)], spacegroup=221, cellpar=[6.03,6.03,6.03,90,90,90])
#pure_tet_23 = crystal(('Cs', 'Pb', 'I', 'I'), basis=[(0.5,0.0,0.5), (0.0,0.0,0.0), (0.0,0.0,0.5), (-0.212,0.288,0.0)], spacegroup=127, cellpar=[8.33,8.33,6.14,90,90,90])
#pure_ortho_23 = crystal(('Cs', 'Pb', 'I', 'I'), basis=[(0.9835,0.4451,0.25), (0.0,0.0,0.5), (0.1973,0.2993,0.5306), (-0.0631,-0.0017,0.75)], spacegroup=62, cellpar=[8.22,8.56,12.01,90,90,90])
#view(pure_cubic)

# 3x3x3, 135-atom cubic supercell
cubic_supercell = make_supercell(pure_cubic_13,[[3,0,0],[0,3,0],[0,0,3]],order="atom-major")
cubic_supercell = [cubic_supercell.repeat((1,1,1))]
write_vasp("cubic_13_333.vasp",cubic_supercell,direct=True)
#view(cubic_supercell)

# 2x2x4, 160-atom tetragonal supercell
tet_supercell = make_supercell(pure_tet_13,[[2,0,0],[0,2,0],[0,0,4]],order="atom-major")
tet_supercell = [tet_supercell.repeat((1,1,1))]
write_vasp("tet_13_224.vasp",tet_supercell,direct=True)

# 2x2x2, 160-atom orthorhombic supercell
ortho_supercell = make_supercell(pure_ortho_13,[[2,0,0],[0,2,0],[0,0,2]],order="atom-major")
ortho_supercell = [ortho_supercell.repeat((1,1,1))]
write_vasp("ortho_13_222.vasp",ortho_supercell,direct=True)
#view(ortho_supercell)

# cluster space: defines the alloy system, max no. sites per cluster, cutoff radii for max cluster diatance
# choose cutoffs such that size of cluster space order of magnitude ~10 (longer-range correlations won't contribute much to total E anyway so don't need to include too many)
# cluster radii output is not max cluster distance but instead average distance of each site to centre of the cluster
# so max radius ~ 1/2 cutoff
cs_cubic_13 = ClusterSpace(pure_cubic_13, [13, 8], [['Cs'], ['Pb'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br']])
cs_tet_13 = ClusterSpace(pure_tet_13, [11, 8], [['Cs'], ['Cs'], ['Pb'], ['Pb'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br']])
cs_ortho_13 = ClusterSpace(pure_ortho_13, [11, 8], [['Cs'], ['Cs'], ['Cs'], ['Cs'], ['Pb'], ['Pb'], ['Pb'], ['Pb'],
['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br'], ['I', 'Br']])

#print(cs_cubic_13)
#print(cs_tet_13)
print(cs_ortho_13)


# target concentrations: should correspond to final SQS compositions A_1-x, B_x
x = 1/3

#generate_sqs(cs_cubic_13, x, cubic_supercell, 100000)
#generate_sqs(cs_tet_13, x, tet_supercell, 100000)
generate_sqs(cs_ortho_13, x, ortho_supercell, 100000)