# -*- coding: utf-8 -*-
"""
@Time : 2024/06/29/ 11:01
@Auth : Scandium, kanglihn@163.com
@File :  xyz2arc.py
@IDE : Pycharm

"""

N_atom = 144
N_step = 25
N_skip = 1
xyzfile_name = "Docker.docker.struc1.allopt.xyz"
arcfile_name = str(xyzfile_name[:-4]) + ".arc"

element = [[] for _ in range(N_step)]
rx = [[] for _ in range(N_step)]
ry = [[] for _ in range(N_step)]
rz = [[] for _ in range(N_step)]

with open(xyzfile_name, "r") as xyzfile:
    lines = xyzfile.readlines()
    for i in range(N_step):
        for j in range(N_atom):
            line_tmp = lines[i * (N_atom + 2) + j + 2].strip()
            line_tmp = line_tmp.split()
            element[i].append(str(line_tmp[0]))
            rx[i].append(float(line_tmp[1]))
            ry[i].append(float(line_tmp[2]))
            rz[i].append(float(line_tmp[3]))

# print(element[1])

# print(rx[0])

with open(arcfile_name, 'w') as arcfile:
    arcfile.write('!BIOSYM archive 3' + '\n')
    arcfile.write('PBC=OFF' + '\n')

    for i in range(N_step):
        if i % N_skip == 0:
            arcfile.write('                       -11571.7839' + '\n')
            arcfile.write("!DATE     Jun 29 10:24:00 2024" + '\n')
            for j in range(N_atom):
                arcfile.write("{:2}".format(element[i][j]) + "   " + "{:15.9f}".format(rx[i][j])
                              + "{:15.9f}".format(ry[i][j]) + "{:15.9f}".format(rz[i][j]) + str(" " * 21)
                              + "{:2}".format(element[i][j]) + "{:7.3f}".format(0.000) + '\n')

            arcfile.write("end" + "\n")
            arcfile.write("end" + "\n")






