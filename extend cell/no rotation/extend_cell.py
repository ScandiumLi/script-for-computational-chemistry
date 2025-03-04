from ase import io, Atoms
import numpy as np

TPA = io.read("TPA.xyz")

# print(TPA.symbols[2], TPA.positions[2])

number_cell = 5
a = 7.730
b = 6.440
c = 3.749
alpha = 92.749 * np.pi / 180
beta = 109.15 * np.pi / 180
gamma = 95.95 * np.pi / 180

vector_a = np.array([a, 0.0, 0.0])
vector_b = np.array([b * np.cos(gamma), b * np.sin(gamma), 0.0])
vector_cx = c * np.cos(beta)
vector_cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
vector_cz = c * np.sqrt((np.sin(beta) * np.sin(gamma))**2 -
                        (np.cos(alpha) - np.cos(beta)*np.cos(gamma))**2)/np.sin(gamma)

vector_c = np.array([vector_cx, vector_cy, vector_cz])

# print (vector_a, vector_b, vector_c)

tmp = []
symb = []
# print(mol_tmp)
for i in range(0, number_cell):
    for j in range(0, number_cell):
        for k in range(0, number_cell):
            for n in range(len(TPA.positions)):
                symb.append(TPA.symbols[n])
                # print(mol_tmp[n])
                trans_vector = (np.array(vector_a) * i + np.array(vector_b) * j
                                + np.array(vector_c) * k)
                tmp.append(TPA.positions[n] + trans_vector)
        # print(tmp[n - 1])

TPA_cluster = Atoms(symb, tmp)
io.write('TPA_cluster.xyz', TPA_cluster, format='xyz')
