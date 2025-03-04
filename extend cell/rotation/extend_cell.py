from ase import io, Atoms
import numpy as np

R_Au3 = io.read("R-Au3-init.xyz")

# print(TPA.symbols[2], TPA.positions[2])

number_cell = 3
a = 28.714/2
b = 28.714/2
c = 28.714/2
alpha = 90 * np.pi / 180
beta = 90 * np.pi / 180
gamma = 90 * np.pi / 180

vector_a = np.array([a, 0.0, 0.0])
vector_b = np.array([b * np.cos(gamma), b * np.sin(gamma), 0.0])
vector_cx = c * np.cos(beta)
vector_cy = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
vector_cz = c * np.sqrt((np.sin(beta) * np.sin(gamma)) ** 2 -
                        (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) ** 2) / np.sin(gamma)

vector_c = np.array([vector_cx, vector_cy, vector_cz])

# print (vector_a, vector_b, vector_c)

initial_point = (R_Au3.positions[0] + R_Au3.positions[51] + R_Au3.positions[102]) / 3

# print(axis_vector)

def rot_matrix(theta, end_point):
    axis_vector = (end_point - initial_point) / (np.linalg.norm(end_point - initial_point))
    rot_mat = [[0.0 for k in range(3)] for j in range(3)]
    ci = np.cos(theta)
    a = 1 - ci
    s = np.sin(theta)
    rot_mat[0][0] = a * axis_vector[0] ** 2 + ci
    rot_mat[0][1] = a * axis_vector[0] * axis_vector[1] - s * axis_vector[2]
    rot_mat[0][2] = a * axis_vector[0] * axis_vector[2] + s * axis_vector[1]
    rot_mat[1][0] = a * axis_vector[0] * axis_vector[1] + s * axis_vector[2]
    rot_mat[1][1] = a * axis_vector[1] ** 2 + ci
    rot_mat[1][2] = a * axis_vector[1] * axis_vector[2] - s * axis_vector[0]
    rot_mat[2][0] = a * axis_vector[0] * axis_vector[2] - s * axis_vector[1]
    rot_mat[2][1] = a * axis_vector[1] * axis_vector[2] + s * axis_vector[0]
    rot_mat[2][2] = a * axis_vector[2] ** 2 + ci

    return rot_mat


tmp = []
symb = []

# print(init_point)
# print(end_point)

# print(mol_tmp)
for i in range(0, number_cell):
    for j in range(0, number_cell):
        for k in range(0, number_cell):
            for n in range(len(R_Au3.positions)):
                symb.append(R_Au3.symbols[n])
                # print(mol_tmp[n])
                end_point = np.array([0.0, 0.0, 0.0])
                R_Au3_new = R_Au3.positions[n]
                trans_vector = (np.array(vector_a) * i + np.array(vector_b) * j
                                + np.array(vector_c) * k)

                end_point_x = initial_point + np.array([10.0, 0.0, 0.0])
                end_point_y = initial_point + np.array([0.0, 10.0, 0.0])
                end_point_z = initial_point + np.array([0.0, 0.0, 10.0])

                R_Au3_new = np.dot(rot_matrix(np.pi * i, end_point_y), (R_Au3_new - initial_point))
                R_Au3_new = np.dot(rot_matrix(np.pi * j, end_point_z), R_Au3_new)
                R_Au3_new = np.dot(rot_matrix(np.pi * k, end_point_x), R_Au3_new)

                tmp.append(list(R_Au3_new + trans_vector))


        # print(tmp[n - 1])
Cl_init = np.array([-7.426, -7.426, -7.426])

for i in range(number_cell+1):
    for j in range(number_cell+1):
        for k in range(number_cell+1):
            trans_vector = (np.array(vector_a) * i + np.array(vector_b) * j
                            + np.array(vector_c) * k)
            symb.append('Cl')
            tmp.append(list(Cl_init + trans_vector))


R_Au3_cluster = Atoms(symb, tmp)
io.write('R_Au3_cluster.xyz', R_Au3_cluster, format='xyz')
