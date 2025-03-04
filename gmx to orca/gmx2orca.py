from os import remove
from os import path
import re
import glob

# Atom section
Atom_ID = []
atom_type = []
element = []
charge = []
eps = []  # kcal/mol
r_min = []  # Ang.
eps_1_4 = []  # kcal/mol
r_min_1_4 = []  # Ang.

sigma = {}
epsilon = {}

# Bonds section
Bond_ID_first = []
Bond_ID_second = []
bond_type = []
bond_length0 = []  # deg.
k_bond = []  # kcal/mol/Ang./Ang.

# Angle section
Angle_ID_first = []
Angle_ID_second = []
Angle_ID_third = []
angle_type = []
theta0 = []  # deg.
k_angle = []  # kcal/mol/rad./rad.

# Dihedral section
Dihedral_ID_first = []
Dihedral_ID_second = []
Dihedral_ID_third = []
Dihedral_ID_fourth = []
dihedral_type = []
delta = []  # deg.
k_dihedral = []  # kcal/mol
multiplicity = []  # n

if path.isfile("ORCAFF.prms"):
    remove("ORCAFF.prms")
with open(glob.glob("TPA.itp")[0], 'r', encoding='utf-8') as txt1:

####################
### match fields ###
####################

    Pattern_atomtype = re.compile(r"\sA\s")
    Pattern_atom = re.compile(r"qtot")
    Pattern_bond = re.compile(r"[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,3}[0-9]{1,2}\.[0-9]{5}\s")
    Pattern_angle = re.compile(r"^\s{2,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,"
                               r"3}[0-9]{1,3}\.[0-9]{7}\s")
    Pattern_dihedral = re.compile(r"^\s{2,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,6}[0-9]{1,5}\s{1,"
                                  r"6}[0-9]{1,5}\s{1,3}[0-9]{1,3}\.[0-9]{7}\s")

##################
### read data ###
##################

    for line in txt1:

        if re.findall(Pattern_atomtype, line):
            atomtype_tmp = line.strip()
            atomtype_tmp = atomtype_tmp.split()

            sigma[atomtype_tmp[0]] = atomtype_tmp[5]
            epsilon[atomtype_tmp[0]] = atomtype_tmp[6]

        # print(line)
        if re.findall(Pattern_atom, line):
            # print(line)
            atom_tmp = line.strip()
            atom_tmp = atom_tmp.split()
            # print(atom_tmp)
            Atom_ID.append(atom_tmp[0])
            atom_type.append(atom_tmp[1])
            if re.findall(r"[A-Za-z]", atom_tmp[4][1:]):
                element.append(atom_tmp[4][:2])
            else:
                element.append(atom_tmp[4][0])
            charge.append(atom_tmp[6])
            eps.append(epsilon[atom_tmp[1]])
            eps_1_4.append(epsilon[atom_tmp[1]])
            r_min.append(sigma[atom_tmp[1]])
            r_min_1_4.append(sigma[atom_tmp[1]])

        if re.findall(Pattern_bond, line):
            # print(line)
            bond_tmp = line.strip()
            bond_tmp = bond_tmp.split()
            Bond_ID_first.append(bond_tmp[0])
            Bond_ID_second.append(bond_tmp[1])
            bond_type.append(bond_tmp[2])
            bond_length0.append(bond_tmp[3])
            k_bond.append(bond_tmp[4])

        if re.findall(Pattern_angle, line):
            # print(line)
            angle_tmp = line.strip()
            angle_tmp = angle_tmp.split()
            Angle_ID_first.append(angle_tmp[0])
            Angle_ID_second.append(angle_tmp[1])
            Angle_ID_third.append(angle_tmp[2])
            angle_type.append(angle_tmp[3])
            theta0.append(angle_tmp[4])
            k_angle.append(angle_tmp[5])

        if re.findall(Pattern_dihedral, line):
            dihedral_tmp = line.strip()
            dihedral_tmp = dihedral_tmp.split()
            Dihedral_ID_first.append(dihedral_tmp[0])
            Dihedral_ID_second.append(dihedral_tmp[1])
            Dihedral_ID_third.append(dihedral_tmp[2])
            Dihedral_ID_fourth.append(dihedral_tmp[3])
            dihedral_type.append(dihedral_tmp[4])
            delta.append(dihedral_tmp[5])
            k_dihedral.append(dihedral_tmp[6])
            multiplicity.append(dihedral_tmp[7])

##################################################################################################
###################################### unit transform ############################################
### In Gromacs, the potential is (1/2)*k*(r-r0)^2, while the potential is k*(r-r0)^2 in ORCA5.###
############ And the unit of k is KJ/mol/nm^2 in Gromacs, but Kcal/mol/A^2 in ORCA5. ############
#################################################################################################

    nm2ang = 10.0
    cal2J = 4.184

    for i in range(len(Atom_ID)):
        eps[i] = float(eps[i]) * (-1) / (cal2J * nm2ang ** 2)
        r_min[i] = float(r_min[i]) * nm2ang
        eps_1_4[i] = float(eps_1_4[i]) * (-1) / (2 * cal2J * nm2ang ** 2)
        r_min_1_4[i] = float(r_min_1_4[i]) * nm2ang

    for i in range(len(Bond_ID_first)):
        bond_length0[i] = float(bond_length0[i]) * nm2ang
        k_bond[i] = float(k_bond[i]) / (2 * cal2J * nm2ang ** 2)

    for i in range(len(Angle_ID_first)):
        k_angle[i] = float(k_angle[i]) / (2 * cal2J)

    for i in range(len(Dihedral_ID_first)):
        k_dihedral[i] = float(k_dihedral[i]) / cal2J

##################
### write data ###
##################

    with open("ORCAFF.prms", 'a') as txt2:
        txt2.write('$fftype\nGromacs\n')
        txt2.write('$atoms\n')
        txt2.write(str(len(Atom_ID)) + ' ' + str(1) + ' ' + str(4) + '\n')
        for i in range(len(Atom_ID)):
            txt2.write("{:>6}".format(Atom_ID[i]) + "{:>4}".format(element[i]) + "{:>14f}".format(float(charge[i]))
                       + "{:>13f}".format(float(eps[i])) + "{:>13f}".format(float(r_min[i]))
                       + "{:>13f}".format(float(eps_1_4[i])) + "{:>13f}".format(float(r_min_1_4[i])) + '\n')

        txt2.write('$bonds\n')
        txt2.write(str(len(Bond_ID_first)) + ' ' + '2' + ' ' + '2' + '\n')
        for i in range(len(Bond_ID_first)):
            txt2.write("{:>6}".format(Bond_ID_first[i]) + "{:>11}".format(Bond_ID_second[i])
                       + "{:>11f}".format(float(bond_length0[i])) + "{:>11f}".format(float(k_bond[i])) + '\n')

        txt2.write('$angles\n')
        txt2.write(str(len(Angle_ID_first)) + ' ' + '3' + ' ' + '2' + '\n')
        for i in range(len(Angle_ID_first)):
            txt2.write("{:>6}".format(Angle_ID_first[i]) + "{:>9}".format(Angle_ID_second[i])
                       + "{:>9}".format(Angle_ID_third[i]) + "{:>13f}".format(float(theta0[i]))
                       + "{:>13f}".format(float(k_angle[i])) + '\n')

        txt2.write('$dihedrals\n')
        txt2.write(str(len(Dihedral_ID_first)) + ' ' + '4' + ' ' + '3' + '\n')
        for i in range(len(Dihedral_ID_first)):
            txt2.write("{:>6}".format(Dihedral_ID_first[i]) + "{:>9}".format(Dihedral_ID_second[i])
                       + "{:>9}".format(Dihedral_ID_third[i]) + "{:>9}".format(Dihedral_ID_fourth[i])
                       + "{:>13f}".format(float(delta[i])) + "{:>13f}".format(float(k_dihedral[i]))
                       + "{:>9}".format(multiplicity[i]) + '\n')

