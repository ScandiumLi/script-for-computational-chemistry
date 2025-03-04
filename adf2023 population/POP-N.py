from os import remove
from os import path

Number = 40  # 绘制图形所需轨道总数 #

# 定义分子中片段所含原子的序号，按需添加 #
frame1 = "1-7,25-30"
frame2 = "13,20,34,37,39,46,47,51,64,70-114,129-139,181-191,197-206,217-218,222,223,227,228,249,250,254,255,262,263,276-280,319-320,333-336,356-359,365,366,382-387,397-407,426,427,448-458,474,475"
frame3 = "8-12,14-19,21-24,31-33,35,36,38,40-45,48,49,50,52-63,65-69,115-128,140-180,192-196,207-216,219,220,221,224-226,229-248,251-253,256-261,264-275,281-318,321-332,337-355,360-364,367-381,388-396,408-425,428-447,459-473,476-489"

# 片段对轨道的贡献 #
sum_1 = []
sum_2 = []
sum_3 = []
sum_total = []


# 用于将分子中原子序号转换为ADF输出文件原子片段序号的函数 #
def n_trans(str2):
    str2_tmp = str2.strip(",")
    str2_tmp = str2_tmp.split(",")
    str1 = []
    n_list = []
    str1.extend(str2_tmp)
    for i in range(len(str1)):
        str1_tmp = str1[i].strip("-")
        str1_tmp = str1_tmp.split("-")
        n_list.extend(list(range(int(str1_tmp[0]), int(str1_tmp[-1]) + 1)))
    new_frame = []
    with open("frame", 'r') as txt1:
        line_tmp = txt1.readlines()
        for i in range(len(line_tmp)):
            line1 = line_tmp[i].strip()
            line1 = line1.split()
            for j in range(len(n_list)):
                if int(line1[2]) == n_list[j]:
                    new_frame.append(line1[0])
                else:
                    continue
    return new_frame


# 序号转换 #
frame1_new = n_trans(frame1)
frame2_new = n_trans(frame2)
frame3_new = n_trans(frame3)

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

# 读取orbits文件内容到number列表中 #
number = []
with open("orbits", 'r') as txt3:
    lines = txt3.readlines()
    num_line = len(lines)
    for i in range(num_line):
        line = lines[i].strip()
        line = line.split()
        number.append(line[-1])

# 初始化sum_1、sum_2、sum_3 #
for i in range(Number):
    sum_1.append(0.0)
    sum_2.append(0.0)
    sum_3.append(0.0)
    sum_total.append(0.0)

# 求和 #
for i in range(Number):
    for j in range(len(pop[0])):
        if number[j] in frame1_new:
            sum_1[i] += pop[i][j]
        elif number[j] in frame2_new:
            sum_2[i] += pop[i][j]
        elif number[j] in frame3_new:
            sum_3[i] += pop[i][j]
        else:
            continue
        sum_total[i] = sum_1[i] + sum_2[i] + sum_3[i]

# 读取Orbital_Energy文件中能量值到energy列表中 #
energy = []
with open("Orbital_Energy", 'r') as txt4:
    lines = txt4.readlines()
    h = 0
    for k in range(len(lines)):
        if '2.00' in lines[k][7:15] and '0.00' in lines[k + 1][7:15]:
            h += k + 1
    for j in range(Number):
        line = lines[h + j - 20].strip()
        line = line.split()
        energy.append(line[0])

# 将结果写入pop.data文件 #
if path.isfile("pop.data"):
    remove("pop.data")
with open("pop.data", 'a') as txt5:
    txt5.write("Energy" + '\t' + "frame1" + '\t' + "frame2" + '\t'  + "frame3" + '\t' + "total" + '\n')
    for i in range(Number):
        txt5.write(energy[i] + '\t' + format(sum_1[i], '-.2f') + '\t' +
                   format(sum_2[i], '-.2f') + '\t' + format(sum_3[i], '-.2f') + '\t' + format(sum_total[i], '-.2f') + '\n')
