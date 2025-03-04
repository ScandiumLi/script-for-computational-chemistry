from ase import io
from ase import Atoms
from ase import Atom
# from ase.visualize import view
# from ase.io.trajectory import Trajectory
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt

Au13_bisPPh3_S0 = io.read("Au13-bisPPh3Cl2-S0-opt.xyz")
Au13_bisAsPh3_S0 = io.read("Au13-bisAsPh3Cl2-S0-opt.xyz")
Au13_bisSbPh3_S0 = io.read("Au13-bisSbPh3Cl2-S0-opt.xyz")

Au13_bisPPh3_S1 = io.read("Au13-bisPPh3Cl2-S1-opt.xyz")
Au13_bisAsPh3_S1 = io.read("Au13-bisAsPh3Cl2-S1-opt.xyz")
Au13_bisSbPh3_S1 = io.read("Au13-bisSbPh3Cl2-S1-opt.xyz")

# view(Au38S20_S0, viewer = "vmd")

bond_Au13_bisPPh3_S0 = {}
bond_Au13_bisAsPh3_S0 = {}
bond_Au13_bisSbPh3_S0 = {}

bond_Au13_bisPPh3_S1 = {}
bond_Au13_bisAsPh3_S1 = {}
bond_Au13_bisSbPh3_S1 = {}

# core_Au13_dege = list(range(1, 12))

core_Au13_dege = list(range(1, 6))
core_Au13_dege.extend(range(88, 93))

comb_Au13_dege = list(combinations(core_Au13_dege, 2))

comb_Au13_center = ((1, 13), (2, 13), (3, 13), (4, 13), (5, 13), (6, 13),(7, 13), (8, 13), (9, 13), (10, 13), (11, 13), (12, 13))

# comb_Au_P = ((1, 7), (2, 8), (3, 9), (4, 10), (5, 11), (6, 13),(88, 94), (89, 95), (90, 96), (91, 97), (92, 98), (93, 99))

comb_Au_P = ((2, 16), (3, 17), (4, 18), (5, 19), (6, 20), (7, 21), (8, 22), (9, 23), (10, 24), (11, 25))

# comb_Au13_center = ((1, 175), (2, 175), (3, 175), (4, 175), (5, 175), (6, 175),(88, 175), (89, 175), (90, 175), (91, 175), (92, 175), (93, 175))

# comb_Au_Cl = ((1, 14), (12, 15))



for i, j in iter(comb_Au_P):
    if Au13_bisPPh3_S0.get_distance(i - 1, j - 1) <= 3.3:
        bond_Au13_bisPPh3_S0[str(i) + '--' + str(j)] = Au13_bisPPh3_S0.get_distance(i - 1, j - 1)
        bond_Au13_bisPPh3_S1[str(i) + '--' + str(j)] = Au13_bisPPh3_S1.get_distance(i - 1, j - 1)
        # bond_Au13_bisAsPh3_S0[str(i) + '--' + str(j)] = Au13_bisAsPh3_S0.get_distance(i - 1, j - 1)
        # bond_Au13_bisSbPh3_S0[str(i) + '--' + str(j)] = Au13_bisSbPh3_S0.get_distance(i - 1, j - 1)


bond_Au13_bisPPh3_S0_new = dict(sorted(bond_Au13_bisPPh3_S0.items(), key=lambda s: s[1]))


#
plt.scatter(bond_Au13_bisPPh3_S0_new.keys(), bond_Au13_bisPPh3_S0_new.values(), marker='o')
plt.scatter(bond_Au13_bisPPh3_S1.keys(), bond_Au13_bisPPh3_S1.values(), marker='*')

# plt.scatter(bond_Au13_bisAsPh3_S0.keys(), bond_Au13_bisAsPh3_S0.values(), marker='*')
# plt.scatter(bond_Au13_bisSbPh3_S0.keys(), bond_Au13_bisSbPh3_S0.values(), marker='*')
plt.ylim(2.0, 3.5)
plt.show()

for item in iter(bond_Au13_bisPPh3_S0_new.keys()):
    print(item, bond_Au13_bisPPh3_S0[item], bond_Au13_bisPPh3_S1[item])
