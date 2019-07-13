import sys
import os
import time
import shutil
import subprocess

SPIN_CENTER = "" #Chemical name of the spin center
INPUT_NAME = "" #Root of the name of the input file
                #(es. Fe2S2, without .inp)
ES_PROGRAM = "orca" #Program for single point calculations, can be
                    #orca of cp2k (not tested)
START_STEP = 1
ROOT = "" #Directory where the calculation is to be run
ORCA_EXE = "" #Path of Orca executable
SPIN_LADDER_EXE = "" #Path of spin_ladder executable
WFS_OPT_PATH = "wfs_opt"
GEO_OPT_PATH = "geo_opt"
SPIN_LADDER_PATH = "spin_ladder"
INTERVAL_TIME_RECHECK = 100
WAITING_TIME_WFS_FINISHED = 100
WAITING_TIME_GEO_FINISHED = 60
SPIN_STATES = [name for name in os.listdir(os.path.join(ROOT, WFS_OPT_PATH,
                                                        str(START_STEP))) if
               os.path.isdir(os.path.join(ROOT, WFS_OPT_PATH,
                                          str(START_STEP), name))]

#Convert ES_PROGRAM to lowercase in order to avoid issues
ES_PROGRAM = ES_PROGRAM.lower()

def parse_energies(iter_step):
    '''Returns list of energies of each spin state in subfolders of
       wfs_opt/iter_step'''
    energies = []
    for i in sorted(os.listdir(os.path.join(ROOT, WFS_OPT_PATH,
                                            str(iter_step)))):
        if os.path.isdir(os.path.join(ROOT, WFS_OPT_PATH,
                                      str(iter_step), i)):
            file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                         str(iter_step), i, INPUT_NAME
                                         + ".out"), "r")
            file_list = file_obj.readlines()
            file_obj.close()
            for j in reversed(file_list):
                if ES_PROGRAM == "cp2k":
                    if "Total energy:" in j:
                        energies.append(j.split()[-1])
                        break
                elif ES_PROGRAM == "orca":
                    if "FINAL SINGLE POINT ENERGY " in j:
                        energies.append(j.split()[-1])
                        break
    return energies

def print_energies(iter_step):
    '''Prints Energies.dat in spin_ladder/iter_step'''
    file_obj = open(os.path.join(ROOT, SPIN_LADDER_PATH, str(iter_step),
                                 "Energies.dat"), "w")
    energies = parse_energies(iter_step)
    for i in energies:
        file_obj.write(i + "\n")
    file_obj.close()

def parse_spin_moments(iter_step):
    '''Returns list of spin moments of spin centers in each spin state
       in subfolders of wfs_opt/iter_step'''
    spin_moments = []
    for i in sorted(os.listdir(os.path.join(ROOT, WFS_OPT_PATH,
                                            str(iter_step)))):
        if os.path.isdir(os.path.join(ROOT, WFS_OPT_PATH,
                                      str(iter_step), i)):
            if ES_PROGRAM == "cp2k":
                file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                             str(iter_step), i,
                                             INPUT_NAME + ".mulliken"),
                                "r")
                file_list = file_obj.readlines()
                file_obj.close()
                spin_moments.append([])
                for j in file_list[5:-3]:
                    if SPIN_CENTER in j:
                        spin_moments[-1].append(j.split()[-1])
                    if "Mulliken Population Analysis" in j:
                        break
            elif ES_PROGRAM == "orca":
                file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                             str(iter_step), i,
                                             INPUT_NAME + ".out"),
                                "r")
                file_list = file_obj.readlines()
                file_obj.close()
                spin_moments.append([])
                for j in reversed(range(len(file_list))):
                    if "Sum of atomic charges" in file_list[j]:
                        bottom_index = j
                    if "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS" \
                        in file_list[j]:
                        top_index = j
                for j in file_list[top_index:bottom_index]:
                    if SPIN_CENTER in j:
                        spin_moments[-1].append(j.split()[-1])
    return spin_moments

def print_spin_moments(iter_step):
    '''Prints M_values.dat in spin_ladder/iter_step'''
    file_obj = open(os.path.join(ROOT, SPIN_LADDER_PATH, str(iter_step),
                                 "M_values.dat"), "w")
    spin_moments = parse_spin_moments(iter_step)
    file_obj.write(str(len(spin_moments)) + " " +
                   str(len(spin_moments[0])) + "\n")
    for i in spin_moments:
        file_obj.write("\t".join(i) + "\n")
    file_obj.close()

def parse_S2(iter_step):
    '''Returns list of S**2 for each spin state in subfolders of
       wfs_opt/iter_step'''
    S2 = []
    for i in sorted(os.listdir(os.path.join(ROOT, WFS_OPT_PATH,
                                            str(iter_step)))):
        if os.path.isdir(os.path.join(ROOT, WFS_OPT_PATH,
                                      str(iter_step), i)):
            file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                         str(iter_step), i, INPUT_NAME
                                         + ".out"), "r")
            file_list = file_obj.readlines()
            file_obj.close()
            for j in reversed(file_list):
                if ES_PROGRAM == "cp2k":
                    if "Ideal and single determinant S**2 :" in j:
                        S2.append(j.split()[-1])
                        break
                elif ES_PROGRAM == "orca":
                    if "Expectation value of <S**2> " in j:
                        S2.append(j.split()[-1])
                        break
    return S2

def print_S2(iter_step):
    '''Prints S2.dat in spin_ladder/iter_step'''
    file_obj = open(os.path.join(ROOT, SPIN_LADDER_PATH,
                                 str(iter_step), "S2_tot.dat"), "w")
    S2 = parse_S2(iter_step)
    for i in S2:
        file_obj.write(i + "\n")
    file_obj.close()

def parse_gradients(iter_step):
    '''Returns list of gradients of each spin state in subfolders of
       wfs_opt/iter_step'''
    gradients = []
    for i in sorted(os.listdir(os.path.join(ROOT, WFS_OPT_PATH,
                                            str(iter_step)))):
        if os.path.isdir(os.path.join(ROOT, WFS_OPT_PATH,
                                      str(iter_step), i)):
            if ES_PROGRAM == "cp2k":
                file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                             str(iter_step), i,
                                             INPUT_NAME + ".xyz"), "r")
                file_list = file_obj.readlines()
                file_obj.close()
                gradients.append([])
                for j in file_list[4:-1]:
                    for k in j.split()[-3:]:
                        gradients[-1].append(str(-1 * float(k)))
            elif ES_PROGRAM == "orca":
                file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                             str(iter_step), i,
                                             INPUT_NAME + ".engrad"),
                                "r")
                file_list = file_obj.readlines()
                file_obj.close()
                gradients.append([])
                for j, k in enumerate(file_list):
                    if "The current gradient in Eh/bohr" in k:
                        top_index = j
                    elif "The atomic numbers and current coordinates \
                        in Bohr" in file_list[j]:
                        bottom_index = j
                for k in file_list[top_index + 2:bottom_index - 1]:
                    gradients[-1].append(k.strip())
    return gradients

def print_gradients(iter_step):
    '''Prints engrad.dat in spin_ladder/iter_step'''
    file_obj = open(os.path.join(ROOT, SPIN_LADDER_PATH, str(iter_step),
                                 "engrad.dat"), "w")
    gradients = parse_gradients(iter_step)
    file_obj.write(str(len(gradients)) + " " + str(len(gradients[0]))
                   + "\n")
    for i in gradients:
        file_obj.write("\t".join(i) + "\n")
    file_obj.close()

def geo_opt_further(iter_step):
    '''Copy INPUT_NAME.xyz from geo_opt to next iter_step of wfs_opt
       but first modify it adding indexes for magnetic atoms'''
    file_obj = open(os.path.join(ROOT, GEO_OPT_PATH, INPUT_NAME
                                 + ".xyz"), "r")
    file_list = file_obj.readlines()
    file_obj.close()
    counter = 1
    for i in range(2, len(file_list)):
        if SPIN_CENTER in file_list[i]:
            file_list[i] = file_list[i].replace(SPIN_CENTER, SPIN_CENTER
                                                + str(counter))
            counter += 1
    for i in SPIN_STATES:
        if ES_PROGRAM == "cp2k":
            file_obj = open(os.path.join(ROOT, WFS_OPT_PATH,
                                         str(iter_step), str(i),
                                         INPUT_NAME + "_last_step.xyz"),
                            "w")
            for j in file_list:
                file_obj.write(j)
            file_obj.close()
        elif ES_PROGRAM == "orca":
            shutil.copy(os.path.join(ROOT, GEO_OPT_PATH, INPUT_NAME
                                     + ".xyz"),
                        os.path.join(ROOT, WFS_OPT_PATH,
                                     str(iter_step), str(i), INPUT_NAME
                                     + "_last_step.xyz"))

def wfs_opt(iter_step):
    '''Run wfs optimisation on all BS states'''
    submit_wfs(iter_step)
    time.sleep(WAITING_TIME_WFS_FINISHED)
    wfs_opt_next(iter_step)
    wfs_opt_further(iter_step)

def submit_wfs(iter_step):
    '''Submits wfs optimisation on all BS states using CP2K'''
    node_names = open(os.path.join(ROOT, 'machinefile'),
                      'r').readlines()
    for i in sorted(SPIN_STATES):
        run_dir = os.path.join(ROOT, WFS_OPT_PATH, str(iter_step),
                               str(i))
        #Print if there is no input file
        if (INPUT_NAME + ".inp") not in os.listdir(run_dir):
            print("Oops! There is not input file for wfs opt")
        else:
            if (i == '0' or i == '6' or i == '5' or i == '7'):
                node_name = node_names[0].rstrip()
                if ES_PROGRAM == "cp2k":
                    subprocess.call("ssh {} 'cd {}; mpirun -np 18 \
                                    cp2k.psmp {}.inp > {}.out 2> \
                                    {}.err &'".format(node_name,
                                                      run_dir,
                                                      INPUT_NAME,
                                                      INPUT_NAME,
                                                      INPUT_NAME),
                                    shell=True)
                elif ES_PROGRAM == "orca":
                    subprocess.call("ssh {} 'cd {}; {} {}.inp -d -v > \
                                    {}.out 2> {}.err'".format(
                                        node_name, run_dir, ORCA_EXE,
                                        INPUT_NAME, INPUT_NAME, INPUT_NAME),
                                    shell=True)
            elif (i == '1' or i == '2' or i == '3' or i == '4'):
                node_name = node_names[1].rstrip()
                if ES_PROGRAM == "cp2k":
                    subprocess.call("ssh {} 'cd {}; mpirun -np 18 \
                                    cp2k.psmp {}.inp > {}.out 2> \
                                    {}.err &'".format(node_name,
                                                      run_dir,
                                                      INPUT_NAME,
                                                      INPUT_NAME,
                                                      INPUT_NAME),
                                    shell=True)
                elif ES_PROGRAM == "orca":
                    subprocess.call("ssh {} 'cd {}; {} {}.inp -d -v > \
                                     {}.out 2> {}.err'".format(
                                         node_name, run_dir, ORCA_EXE,
                                         INPUT_NAME, INPUT_NAME, INPUT_NAME),
                                    shell=True)

def wfs_opt_next(iter_step):
    '''Creates folder into next iter_step, copy input files and optimized MOs'''
    for i in SPIN_STATES:
        while not wfs_finished_succeed(os.path.join(ROOT, WFS_OPT_PATH,
                                                    str(iter_step),
                                                    str(i), INPUT_NAME
                                                    + ".out")):
            time.sleep(INTERVAL_TIME_RECHECK)
        if os.path.isdir(os.path.join(ROOT, WFS_OPT_PATH,
                                      str(iter_step + 1), str(i))):
            shutil.rmtree(os.path.join(ROOT, WFS_OPT_PATH,
                                       str(iter_step + 1), str(i)))
        os.makedirs(os.path.join(ROOT, WFS_OPT_PATH, str(iter_step + 1),
                                 str(i)))
        if ES_PROGRAM == "cp2k":
            shutil.copy(os.path.join(ROOT, WFS_OPT_PATH, str(iter_step),
                                     str(i), INPUT_NAME + ".wfn"),
                        os.path.join(ROOT, WFS_OPT_PATH, str(iter_step
                                                             + 1),
                                     str(i), INPUT_NAME
                                     + "_last_step.wfn"))
        elif ES_PROGRAM == "orca":
            shutil.copy(os.path.join(ROOT, WFS_OPT_PATH, str(iter_step),
                                     str(i), INPUT_NAME + ".gbw"),
                        os.path.join(ROOT, WFS_OPT_PATH, str(iter_step
                                                             + 1),
                                     str(i), INPUT_NAME
                                     + "_last_step.gbw"))
        shutil.copy(os.path.join(ROOT, WFS_OPT_PATH, str(iter_step),
                                 str(i), INPUT_NAME + ".inp"),
                    os.path.join(ROOT, WFS_OPT_PATH, str(iter_step + 1), str(i),
                                 INPUT_NAME + ".inp"))

def wfs_opt_further(iter_step):
    '''Generates M_values.dat, S2.dat, Energies.dat and engrad.dat in
       spin_ladder/iter_step'''
    if (os.path.isdir(os.path.join(ROOT, SPIN_LADDER_PATH,
                                   str(iter_step)))):
        shutil.rmtree(os.path.join(ROOT, SPIN_LADDER_PATH,
                                   str(iter_step)))
    os.makedirs(os.path.join(ROOT, SPIN_LADDER_PATH, str(iter_step)))
    print_energies(iter_step)
    print_spin_moments(iter_step)
    print_S2(iter_step)
    print_gradients(iter_step)

def wfs_finished_succeed(fullname):
    '''Return True if CP2K wfs optimization has finished'''
    file_obj = open(fullname, "r")
    file_list = file_obj.readlines()
    file_obj.close()
    if ES_PROGRAM == "cp2k":
        return WFS_OPT_PATH in file_list[-1] #WFS_OPT_PATH is because CP2K
                                             #prints the path of the
                                             #running directory in his last
                                             #line
    elif ES_PROGRAM == "orca":
        return "ORCA TERMINATED NORMALLY" in file_list[-2]
    return False

def spin_ladder(iter_step):
    '''Executes SPIN_LADDER_EXE'''
    spin_ladder_dir = os.path.join(ROOT, SPIN_LADDER_PATH,
                                   str(iter_step))
    os.chdir(spin_ladder_dir)
    spin_ladder_run = subprocess.call("{} < {}/heisenberg.inp > {} 2> \
                                      {}".format(SPIN_LADDER_EXE, ROOT,
                                                 os.path.join(
                                                     spin_ladder_dir,
                                                     "spin_ladder.out"),
                                                 os.path.join(
                                                     spin_ladder_dir,
                                                     "spin_ladder.err")),
                                      shell=True)
    if spin_ladder_run != 0:
        raise Exception(spin_ladder_run)
    else:
        os.chdir(ROOT)
        spin_ladder_further(iter_step)

def spin_ladder_further(iter_step):
    '''Copies GS.extcomp.out to geo_opt'''
    while not os.path.isfile(os.path.join(ROOT, SPIN_LADDER_PATH,
                                          str(iter_step),
                                          "GS.extcomp.out")):
        time.sleep(INTERVAL_TIME_RECHECK)
    shutil.copy(os.path.join(ROOT, SPIN_LADDER_PATH,
                             str(iter_step), "GS.extcomp.out"),
                os.path.join(ROOT, GEO_OPT_PATH, INPUT_NAME
                             + ".extcomp.out"))

def main():
    '''Main optimisation step'''
    iter_step_file = open(os.path.join(ROOT, GEO_OPT_PATH, 'iter_step'),
                          'r')
    iter_step = int(iter_step_file.readline())
    iter_step_file.close()
    geo_opt_further(iter_step)
    print("Preparing geometry for step " + str(iter_step) + "\n")
    wfs_opt(iter_step)
    print("HS and BS wfs optmization of step " + str(iter_step) + "\n")
    spin_ladder(iter_step)
    print("Spin Ladder of step               " + str(iter_step) + "\n")
    sys.stdout.flush()
    iter_step += 1
    iter_step_file = open(os.path.join(ROOT, GEO_OPT_PATH, 'iter_step'),
                          'w')
    iter_step_file.write(str(iter_step))
    iter_step_file.close()

if __name__ == "__main__":
    main()
