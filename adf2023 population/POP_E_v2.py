from os import remove
from os import path
import re
import glob
import atexit


Number = 40  # 绘制图形所需轨道总数 #
# 输入金属元素
metal_element = ["Au"]
# 输入非金属元素, 不包含CH
nonmetal_element = ["N", "Br"]
# ADF计算结果文件
filename = "pop.out"


N_start = []
if path.isfile("orbits"):
    remove("orbits")
with open(glob.glob("pop.out")[0], 'r', encoding='utf-8') as txt1:
    with open("orbits", 'a') as txt2:
        Pattern2 = re.compile(r"[0-9]{3}\sau\s\s[A-Za-z]{1,2}")
        for line in txt1:
            # print(line)
            if re.findall(Pattern2, line):
                line_tmp = line.strip()
                line_tmp = line_tmp.split()
                # print(line_tmp)
                N_start.append(line_tmp[1])
                txt2.write(str(line_tmp[5]) + '\t' + str(line_tmp[8][0:1]) + '\t' + str(line_tmp[-1]) + ' \n')
            else:
                continue

if path.isfile("populations"):
    remove("populations")
with open(glob.glob("pop.out")[0], 'r',encoding='utf-8') as txt1:
    with open("populations", 'a') as txt2:
        lines = txt1.readlines()
        for i in range(len(lines)):
            if "SFO contributions (%) per orbital" in lines[i]:
                n = 0
                for k in range(len(lines) - i):
                    if (':' in lines[i + k + 7]):
                        n += 1
                    else:
                        break
                for j in range(n):
                    line1 = lines[i + j + 7].strip(":")
                    line1 = line1.split(":")
                    # print(N_start[0])
                    if int(line1[0]) >= int(N_start[0]):
                        # txt2.write(lines[i + j + 7][10:-1] + '\t')
                        p = Number//14
                        # print(p)
                        for t in range(p+1):
                            txt2.write(lines[i + j + t * (n + 5) + 7][10:-1] + '\t')
                        txt2.write('\n')
                    else:
                        continue

if (path.isfile("Orbital_Energy")):
    remove("Orbital_Energy")
with open(glob.glob("pop.out")[0], 'r', encoding='utf-8') as txt1:
    with open("Orbital_Energy", 'a') as txt2:
        # lines = txt1.readlines()
        Pattern1 = re.compile(r"[0-9]{2}\s{2,6}[0-9]{1,3}\sA\s")
        for line in txt1:
            # print(line)
            if re.findall(Pattern1, line):
                line_tmp = line.strip()
                line_tmp = line_tmp.split()
                txt2.write(str(line_tmp[0]) + '\t\t' + str(line_tmp[1]) + '\t' + str(line_tmp[2]) + ' ' + str(
                    line_tmp[3]) + ' \n')
            else:
                continue



def generate_metal_lists(element_list, init_length):
    """
    生成包含初始化浮点数的列表
    :param metal_list: 金属列表，如 ["Au", "Ag"]
    :param init_length: 初始化列表长度（每个列表填充 init_length 个 0.0）
    :return: 包含生成列表的字典，同时动态注入全局变量
    """
    generated = {}
    for element in element_list:
        # 创建两个初始化列表
        sp_list = [0.0] * init_length
        d_list = [0.0] * init_length
        
        # 存储到字典
        generated[f"sum_{element}_SP"] = sp_list
        generated[f"sum_{element}_D"] = d_list
        
        # 动态注入全局变量（可选）
        globals()[f"sum_{element}_SP"] = sp_list
        globals()[f"sum_{element}_D"] = d_list
    return generated


def generate_nonmetal_lists(element_list, init_length):
    """
    生成包含初始化浮点数的列表
    :param metal_list: 金属列表，如 ["Au", "Ag"]
    :param init_length: 初始化列表长度（每个列表填充 init_length 个 0.0）
    :return: 包含生成列表的字典，同时动态注入全局变量
    """
    generated = {}
    for element in element_list:
        # 创建两个初始化列表
        elem_list = [0.0] * init_length
        
        # 存储到字典
        generated[f"sum_{element}"] = elem_list
       
        # 动态注入全局变量（可选）
        globals()[f"sum_{element}"] = elem_list
    return generated


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

# 初始化#
metal_dict = generate_metal_lists(metal_element, Number)
nonmetal_dict = generate_nonmetal_lists(nonmetal_element, Number)
sum_CH = [0.0] * Number
sum_total = [0.0] * Number


# 求和 #
for i in range(Number):
    for j in range(len(pop[0])):
        for k in range(len(metal_element)):
            if element[j] == metal_element[k] and basis_function[j] == 'S':
                metal_dict[f"sum_{metal_element[k]}_SP"][i] += pop[i][j]
                sum_total[i] += pop[i][j]
            elif element[j] == metal_element[k] and basis_function[j] == 'P':
                metal_dict[f"sum_{metal_element[k]}_SP"][i] += pop[i][j]
                sum_total[i] += pop[i][j] 
            elif element[j] == metal_element[k] and basis_function[j] == 'D':
                metal_dict[f"sum_{metal_element[k]}_D"][i] += pop[i][j]
                sum_total[i] += pop[i][j]
            else:
                continue
        for k in range(len(nonmetal_element)):
            if element[j] == nonmetal_element[k]:
                nonmetal_dict[f"sum_{nonmetal_element[k]}"][i] += pop[i][j]
                sum_total[i] += pop[i][j]
            else:
                continue
        if element[j] == 'C':
            sum_CH[i] += pop[i][j]
            sum_total[i] += pop[i][j]
        elif element[j] == 'H':
            sum_CH[i] += pop[i][j]
            sum_total[i] += pop[i][j]
        else:
            continue

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
    title = "Energy"
    for elem in metal_element:
        title += '\t' + f"sum_{elem}_SP" + '\t' + f"sum_{elem}_D"
    for elem in nonmetal_element:
        title += '\t' + f"sum_{elem}" 
    txt5.write(title + '\t' + "sum_CH" + '\t' + "total" + '\n')
    for i in range(Number):
        data = str(energy[i])
        for elem in metal_element:
            data += '\t' + format(metal_dict[f"sum_{elem}_SP"][i], '-.2f') + '\t' + format(metal_dict[f"sum_{elem}_D"][i], '-.2f')
        for elem in nonmetal_element:
            data += '\t' + format(nonmetal_dict[f"sum_{elem}"][i], '-.2f') 
        txt5.write(data + '\t' + format(sum_CH[i], '-.2f') + '\t' + format(sum_total[i], '-.2f') + '\n')

def cleanup():
    remove("frame")
    remove("Orbital_Energy")
    remove("orbits")
    remove("populations")

if __name__ == "__main__":
    atexit.register(cleanup)