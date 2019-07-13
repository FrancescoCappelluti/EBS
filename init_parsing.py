'''Parse BS states outfiles'''
#Main output file must be called INPUT_NAME.out
#Mulliken population file must be called INPUT_NAME.mulliken
#Gradient file must be generated with FORMAT XYZ and be called
#INPUT_NAME.xyz

import os
SPIN_CENTER = "" #Chemical name of the spin center
INPUT_NAME = "" #Root of the name of the input file
                #(es. Fe2S2, without .inp)
ES_PROGRAM = "" #Program for single point calculations, can be
                    #orca of cp2k (not tested)

def parse_energies():
    '''Returns list of energies of each spin state in subfolders of
    wfs_opt/iter_step'''
    energies = []
    for i in sorted(os.listdir(".")):
        if os.path.isdir(os.path.join(".", i)):
            file_obj = open(os.path.join(".", i, INPUT_NAME
                                         + ".out"), "r")
            file_list = file_obj.readlines()
            file_obj.close()
            for j in reversed(file_list):
                if ES_PROGRAM == "orca":
                    if "FINAL SINGLE POINT ENERGY " in j:
                        energies.append(j.split()[-1])
                        break
                elif ES_PROGRAM == "cp2k":
                    if "Total energy:" in j:
                        energies.append(j.split()[-1])
                        break
    return energies

def print_energies():
    '''Print BS states energies in Energies.dat'''
    file_obj = open("Energies.dat", "w")
    energies = parse_energies()
    for i in energies:
        file_obj.write(i + "\n")
    file_obj.close()

def parse_spin_moments():
    '''Returns list of spin moments of spin centers in each spin state
    in subfolders of wfs_opt/iter_step'''
    spin_moments = []
    for i in sorted(os.listdir(".")):
        if os.path.isdir(os.path.join(".", i)):
            if ES_PROGRAM == "orca":
                file_obj = open(os.path.join(".", i, INPUT_NAME
                                             + ".out"), "r")
                file_list = file_obj.readlines()
                file_obj.close()
                spin_moments.append([])
                for j in reversed(range(len(file_list))):
                    if "Sum of atomic charges" in file_list[j]:
                        bottom_index = j
                    if "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS" \
                        in file_list[j]:
                        top_index = j
                        break
                for j in file_list[top_index:bottom_index]:
                    if SPIN_CENTER in j:
                        spin_moments[-1].append(j.split()[-1])
            elif ES_PROGRAM == "cp2k":
                file_obj = open(os.path.join(".", i, INPUT_NAME
                                             + ".mulliken"), "r")
                file_list = file_obj.readlines()
                file_obj.close()
                spin_moments.append([])
                for j in file_list[5:-3]:
                    if SPIN_CENTER in j:
                        spin_moments[-1].append(j.split()[-1])
                    if "Mulliken Population Analysis" in j:
                        break
    return spin_moments

def print_spin_moments():
    '''Print BS states spin moments in M_values.dat'''
    file_obj = open("M_values.dat", "w")
    spin_moments = parse_spin_moments()
    file_obj.write(str(len(spin_moments)) + " "
                   + str(len(spin_moments[0])) + "\n")
    for i in spin_moments:
        file_obj.write("\t".join(i) + "\n")
    file_obj.close()

def parse_S2():
    '''Returns list of S**2 for each spin state in subfolders of
    wfs_opt/iter_step'''
    S2 = []
    for i in sorted(os.listdir(".")):
        if os.path.isdir(os.path.join(".", i)):
            file_obj = open(os.path.join(".", i, INPUT_NAME
                                         + ".out"), "r")
            file_list = file_obj.readlines()
            file_obj.close()
            for j in reversed(file_list):
                if ES_PROGRAM == "orca":
                    if "Expectation value of <S**2> " in j:
                        S2.append(j.split()[-1])
                        break
                elif ES_PROGRAM == "cp2k":
                    if "Ideal and single determinant S**2 :" in j:
                        S2.append(j.split()[-1])
                        break
    return S2

def print_S2():
    '''Print BS states S**2 in S2_tot.dat'''
    file_obj = open("S2_tot.dat", "w")
    S2 = parse_S2()
    for i in S2:
        file_obj.write(i + "\n")
    file_obj.close()

def parse_gradients():
    '''Returns list of gradients of each spin state in subfolders of
    wfs_opt/iter_step'''
    gradients = []
    for i in sorted(os.listdir(".")):
        if os.path.isdir(os.path.join(".", i)):
            if ES_PROGRAM == "orca":
                file_obj = open(os.path.join(".", i, INPUT_NAME +
                                             ".engrad"), "r")
                file_list = file_obj.readlines()
                file_obj.close()
                gradients.append([])
                for j, k in enumerate(file_list):
                    if "The current gradient in Eh/bohr" in k:
                        top_index = j
                    elif "The atomic numbers and current coordinates" \
                        " in Bohr" in k:
                        bottom_index = j
                for k in file_list[top_index + 2:bottom_index - 1]:
                    gradients[-1].append(k.strip())
            elif ES_PROGRAM == "cp2k":
                file_obj = open(os.path.join(".", i, INPUT_NAME
                                             + ".xyz"), "r")
                file_list = file_obj.readlines()
                file_obj.close()
                gradients.append([])
                for j in file_list[4:-1]:
                    for k in j.split()[-3:]:
                        gradients[-1].append(str(-1 * float(k)))
    return gradients

def print_gradients():
    '''Print BS states gradients in engrad.dat'''
    file_obj = open("engrad.dat", "w")
    gradients = parse_gradients()
    file_obj.write(str(len(gradients)) + " " + str(len(gradients[0]))
                   + "\n")
    for i in gradients:
        file_obj.write("\t".join(i) + "\n")
    file_obj.close()

def main():
    '''Parse BS states output files'''
    print_energies()
    print_spin_moments()
    print_S2()
    print_gradients()

if __name__ == "__main__":
    main()
