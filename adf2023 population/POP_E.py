from os import remove
from os import path

Number = 40  # 绘制图形所需轨道总数 #

# 基函数对轨道的贡献 #
sum_Au_SP = []
sum_Au_D = []
sum_S = []
sum_CH = []
sum_total = []

# 读取populations文件内容到pop矩阵 #
with open("populations", 'r') as txt2:
    lines = txt2.readlines()
    num_line = len(lines)
    pop = [[0 for i in range(num_line)] for j in range(Number)]
    for i in range(num_line):
        line = lines[i].strip()
        line = line.split()
        for j in range(Number):
            pop[j][i] += float(line[j])

# 读取orbits文件内容到basis_function列表中 #
basis_function = []
with open("orbits", 'r') as txt3:
    lines = txt3.readlines()
    num_line = len(lines)
    for i in range(num_line):
        line = lines[i].strip()
        line = line.split()
        basis_function.append(line[-2])

element = []
with open("orbits", 'r') as txt3:
    lines = txt3.readlines()
    num_line = len(lines)
    for i in range(num_line):
        line = lines[i].strip()
        line = line.split()
        element.append(line[0])

# 初始化sum_1、sum_2 #
for i in range(Number):
    sum_Au_SP.append(0.0)
    sum_Au_D.append(0.0)
    sum_S.append(0.0)
    sum_CH.append(0.0)
    sum_total.append(0.0)

# 求和 #
for i in range(Number):
    for j in range(len(pop[0])):
        if element[j] == 'Au' and basis_function[j] == 'S':
            sum_Au_SP[i] += pop[i][j]
        elif element[j] == 'Au' and basis_function[j] == 'P':
            sum_Au_SP[i] += pop[i][j]
        elif element[j] == 'Au' and basis_function[j] == 'D':
            sum_Au_D[i] += pop[i][j]
        elif element[j] == 'S':
            sum_S[i] += pop[i][j]
        elif element[j] == 'C':
            sum_CH[i] += pop[i][j]
        elif element[j] == 'H':
            sum_CH[i] += pop[i][j]
        else:
            continue
        sum_total[i] = sum_Au_SP[i] + sum_Au_D[i] + sum_S[i] + sum_CH[i]

# 读取Orbital_Energy文件中能量值到energy列表中 #
energy = []
with open("Orbital_Energy", 'r') as txt4:
    lines = txt4.readlines()
    h = 0
    for k in range(len(lines)):
        if '2.00' in lines[k][7:15] and '0.00' in lines[k + 1][7:15]:
            h += k + 1
    for j in range(Number):
        line = lines[h + j - int(Number/2)].strip()
        line = line.split()
        energy.append(line[0])

# 将结果写入pop.data文件 #
if path.isfile("pop.data_2"):
    remove("pop.data_2")
with open("pop.data_2", 'a') as txt5:
    txt5.write("Energy" + '\t' + "Au_SP" + '\t' + "Au_D" + '\t' +
               "S" + '\t' + "CH" + '\t' + "total" + '\n')
    for i in range(Number):
        txt5.write(energy[i] + '\t' + format(sum_Au_SP[i], '-.2f') + '\t' + format(sum_Au_D[i], '-.2f') + '\t' +
                   format(sum_S[i], '-.2f') + '\t' +
                   format(sum_CH[i], '-.2f') + '\t' + format(sum_total[i], '-.2f') + '\n')
