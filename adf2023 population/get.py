from os import remove
from os import path
import re
import glob

Number = 40  # 绘制图形所需轨道总数,只需更改它 #

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
with open(glob.glob("*.out")[0], 'r', encoding='utf-8') as txt1:
    with open("Orbital_Energy", 'a') as txt2:
        # lines = txt1.readlines()
        Pattern1 = re.compile(r"[0-9]{2}\s{2,6}[0-9]{1,3}\sA\s")
        for line in txt1:
            print(line)
            if re.findall(Pattern1, line):
                line_tmp = line.strip()
                line_tmp = line_tmp.split()
                txt2.write(str(line_tmp[0]) + '\t\t' + str(line_tmp[1]) + '\t' + str(line_tmp[2]) + ' ' + str(
                    line_tmp[3]) + ' \n')
            else:
                continue

if path.isfile("frame"):
    remove("frame")
with open(glob.glob("*.out")[0], 'r', encoding='utf-8') as txt1:
    with open("frame", 'a') as txt2:
        Pattern2 = re.compile(r"[0-9]{1,3}\s\s[A-Za-z]{1,2}\s{46,51}[0-9]{1,3}\s\s[A-Za-z]{1,2}")
        for line in txt1:
            # print(line)
            if re.findall(Pattern2, line):
                line_tmp = line.strip()
                line_tmp = line_tmp.split()
                print(line_tmp)
                txt2.write(str(line_tmp[0]) + '\t' + str(line_tmp[1]) + '\t'
                           + str(line_tmp[2]) + '\t' + str(line_tmp[3]) + ' \n')
            else:
                continue