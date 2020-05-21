#Main output file must be called input_name.out
#Mulliken population file must be called input_name.mulliken
#Gradient file must be generated with FORMAT XYZ and be called input_name.xyz

import os
#To be set
spin_center = ""
input_name = ""
es_program = ""
calc_type = ""

def parse_energies():
        #Returns list of energies of each spin state in subfolders of wfs_opt/iter_step
        energies = []
        for i in sorted(os.listdir(".")):
                if os.path.isdir(os.path.join(".", i)):
                        file_obj = open(os.path.join(".", i, input_name + ".out"), "r")
                        file_list = file_obj.readlines()
                        file_obj.close()
                        for j in reversed(file_list):
                                if es_program == "cp2k":
                                        if "Total energy:" in j:
                                                energies.append(j.split()[-1])
                                                break
                                elif es_program == "orca":
                                        if "FINAL SINGLE POINT ENERGY " in j:
                                                energies.append(j.split()[-1])
                                                break
        return energies

def print_energies():
	file_obj = open("Energies.dat", "w")
	energies = parse_energies()
	for i in energies:
		file_obj.write(i + "\n")
	file_obj.close()

def parse_spin_moments():
        #Returns list of spin moments of spin centers in each spin state in subfolders of wfs_opt/iter_step
        spin_moments = []
        for i in sorted(os.listdir(".")):
                if os.path.isdir(os.path.join(".", i)):
                        if es_program == "cp2k":
                                file_obj = open(os.path.join(".", i, input_name + ".mulliken"), "r")
                                file_list = file_obj.readlines()
                                file_obj.close()
                                spin_moments.append([])
                                for j in file_list[5:-3]:
                                        if spin_center in j:
                                                spin_moments[-1].append(j.split()[-1])
                                        if "Mulliken Population Analysis" in j:
                                                break
                        elif es_program == "orca":
                                file_obj = open(os.path.join(".", i, input_name + ".out"), "r")
                                file_list = file_obj.readlines()
                                file_obj.close()
                                spin_moments.append([])
                                if calc_type == "MP2":
                                    vpattern = "MULLIKEN ATOMIC CHARGES AND SPIN DENSITIES"
                                elif calc_type == "DFT":
                                    vpattern = "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS"
                                for j in reversed(range(len(file_list))):
                                        if "Sum of atomic charges" in file_list[j]:
                                                bottom_index = j
                                        if vpattern in file_list[j]:
                                                top_index = j
                                                break
                                for j in file_list[top_index:bottom_index]:
                                        if spin_center in j:
                                                spin_moments[-1].append(j.split()[-1])
        return spin_moments

def print_spin_moments():
	file_obj = open("M_values.dat", "w")
	spin_moments = parse_spin_moments()
	file_obj.write(str(len(spin_moments)) + " " + str(len(spin_moments[0])) + "\n")
	for i in spin_moments:
		file_obj.write("\t".join(i) + "\n")
	file_obj.close()

def parse_S2():
        #Returns list of S**2 for each spin state in subfolders of wfs_opt/iter_step
        S2 = []
        for i in sorted(os.listdir(".")):
                if os.path.isdir(os.path.join(".", i)):
                        file_obj = open(os.path.join(".", i, input_name + ".out"), "r")
                        file_list = file_obj.readlines()
                        file_obj.close()
                        for j in reversed(file_list):
                                if es_program == "cp2k":
                                        if "Ideal and single determinant S**2 :" in j:
                                                S2.append(j.split()[-1])
                                                break
                                elif es_program == "orca":
                                        if "Expectation value of <S**2> " in j:
                                                S2.append(j.split()[-1])
                                                break
        return S2

def print_S2():
	file_obj = open("S2_tot.dat", "w")
	S2 = parse_S2()
	for i in S2:
		file_obj.write(i + "\n")
	file_obj.close()

def parse_gradients():
        #Returns list of gradients of each spin state in subfolders of wfs_opt/iter_step
        gradients = []
        for i in sorted(os.listdir(".")):
                if os.path.isdir(os.path.join(".", i)):
                        if es_program == "cp2k":
                                file_obj = open(os.path.join(".", i, input_name + ".xyz"), "r")
                                file_list = file_obj.readlines()
                                file_obj.close()
                                gradients.append([])
                                for j in file_list[4:-1]:
                                        for k in j.split()[-3:]:
                                                gradients[-1].append(str(-1 * float(k)))
                        elif es_program == "orca":
                                file_obj = open(os.path.join(".", i, input_name + ".engrad"), "r")
                                file_list = file_obj.readlines()
                                file_obj.close()
                                gradients.append([])
                                for j in range(len(file_list)):
                                        if "The current gradient in Eh/bohr" in file_list[j]:
                                                top_index = j
                                        elif "The atomic numbers and current coordinates in Bohr" in file_list[j]:
                                                bottom_index = j
                                for k in file_list[top_index + 2:bottom_index - 1]:
                                                gradients[-1].append(k.strip())
        return gradients

def print_gradients():
	file_obj = open("engrad.dat", "w")
	gradients = parse_gradients()
	file_obj.write(str(len(gradients)) + " " + str(len(gradients[0])) + "\n")
	for i in gradients:
		file_obj.write("\t".join(i) + "\n")
	file_obj.close()

def main():
	print_energies()
	print_spin_moments()
	print_S2()
	print_gradients()

if __name__ == "__main__":
	main()
