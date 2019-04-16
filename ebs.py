#!/usr/bin/env python
'''

Orca call orca_External to run python script.

hierarchy directories  to run are as following:

* spin_ladder
* geo_opt
** step_init # do normal geometry to generate xyz and gbw file for following wfs_opt, do this by hand
** 0 #optmization iter_step 0 of main loop
** 1 #opt iter_step 1 of main loop

* wfs_opt #only  BS and HS states are permited in this path not any other directories.
** 0 #optmization iter_step 0 of main loop
*** 0 #HS spin state wfs opt
*** 1 #BS1 spin state wfs opt
...
** 1 #opt iter_step 1 of main loop
...


Precedures:
* orca_External
0. get the converged BS states
1. start_step and total steps
2. wfs_opt_path and input files including inp.inp, inp_last_step.xyz and inp_last_step.gbw located in wfs_opt_path/start_step.
   It is better to delete output files from last step.
   Programming started at wfs_opt_path.
4. geo_opt_path and input files including inp.inp and inp_last_step.xyz located in geo_opt_path/start_step
5. out_file
6. inp_file
7. interval_time_recheck
8. waiting_time_wfs_finished
9. cp orca_wfs.sh spin_ladder.sh and heisenberg.inp to $ROOT


* HS state geo optmization under geo_opt/step_init/, cp gbw and xyz to wfs_opt/0/0/...wfs_opt/0/7/
* BS and HS states inputs under such as wfs_opt/0/0/.../wfs_opt/0/7/ for Fe4S4 case


* heisenberg.inp
set the spin moment
* orca_S
* orca_External
* input file and geometry for geo_opt and wfs_opt
* iter_step

'''
import os
import sys
import os.path as op
import delegator
import shutil
import time
import fnmatch
import argparse
from pyparsing import Word, alphas, nums, Regex, Literal, OneOrMore, Combine, CaselessLiteral, Optional, Forward, Or, LineEnd, SkipTo

start_step = 1 
root ="maindirectory"
ORCA_EXE = "pathtoorcaexecutable"
SPIN_LADDER_EXE="pathtospin_ladder.x"
wfs_opt_path = "wfs_opt"
geo_opt_path = "geo_opt"
spin_ladder_path = "spin_ladder"
out_file = "Fe4S4.out"
inp_file = "Fe4S4.inp"
err_file = "out.err"
clean_output = True  # clean output files in wfs_opt
orca_finished_keywords = "ORCA TERMINATED NORMALLY"
interval_time_recheck = 100
waiting_time_wfs_finished = 100
waiting_time_geo_finished = 60
bs_spin_state = 8
spin_state_HS = 0
spin_states = [
    name
    for name in os.listdir(op.join(root, wfs_opt_path, str(start_step)))
    if op.isdir(op.join(root, wfs_opt_path, str(start_step), name))
]

#########################################
#parsing section
SpinCenter = 'Fe'
# you can set the path of files manually
orca_out_files_custom = []  # default under $root
orca_engrad_files_custom = []

##########################################


def main():
    #

    iter_step_file = open(op.join(root, geo_opt_path, 'iter_step'), 'r')
    iter_step = int(iter_step_file.readline())

    geo_opt_further(iter_step)
    print("Preparing geometry for step " + str(iter_step) + "\n")
    wfs_opt(iter_step)
    print("HS and BS wfs optmization of step " + str(iter_step) + "\n")
    spin_ladder(iter_step)
    print("Spin Ladder of step               " + str(iter_step) + "\n")
    sys.stdout.flush()

    iter_step = iter_step + 1
    iter_step_file = open(op.join(root, geo_opt_path, 'iter_step'), 'w')
    iter_step_file.write("%s" % iter_step)
    iter_step_file.close()


def geo_opt_further(iter_step):
    # copy xyz from geo_opt to next iter_step of wfs_opt
    for i in spin_states:
        shutil.copy(
            op.join(root, geo_opt_path, inp_file.split('.')[0] + ".xyz"),
            op.join(root, wfs_opt_path, str(iter_step), str(i),
                    inp_file.split('.')[0] + "_last_step.xyz"))


def spin_ladder(iter_step):
    #####################################################
    spin_ladder_dir = op.join(root, spin_ladder_path, str(iter_step))
    spin_ladder_run = delegator.run('cd ' + spin_ladder_dir  + ';' + SPIN_LADDER_EXE + ' <  ' + root + '/heisenberg.inp > spin_ladder.out 2> spin_ladder.err') 
    if (spin_ladder_run.return_code != 0):
        raise Exception(spin_ladder_run.std_err)
        sys.stdout.flush()
    else:
        spin_ladder_further(iter_step)


        # spin_ladder_next has been done by wfs_opt_further function.
def spin_ladder_next(iter_step):
    if os.path.isdir(op.join(root, spin_ladder_path, str(iter_step + 1))):
        shutil.rmtree(op.join(root, spin_ladder_path, str(iter_step + 1)))
    os.makedirs(op.join(root, spin_ladder_path, str(iter_step + 1)))


def spin_ladder_further(iter_step):
    while not os.path.isfile(op.join(root, spin_ladder_path, str(iter_step),
                                     "GS.extcomp.out")):
        time.sleep(interval_time_recheck)
    else:
        shutil.copy(
            op.join(root, spin_ladder_path, str(iter_step), "GS.extcomp.out"),
            op.join(root, geo_opt_path,
                    inp_file.split('.')[0] + ".extcomp.out"))


def wfs_opt(iter_step):
    #####################################################
    # wfs opt
    #run wfs opt on all BS and HS states
    #clean output files from last steps
    submit_orca_wfs(iter_step)
    time.sleep(waiting_time_wfs_finished)
    wfs_opt_next(iter_step)
    wfs_opt_further(iter_step)


def wfs_opt_next(iter_step):
    for i in spin_states:
        while not orca_finished_succeed(op.join(root, wfs_opt_path, str(
                iter_step), str(i), out_file)):
            time.sleep(interval_time_recheck)

        #create folder into next iter_step and copy input files
        if os.path.isdir(op.join(root, wfs_opt_path, str(iter_step + 1), str(
                i))):
            shutil.rmtree(op.join(root, wfs_opt_path, str(iter_step + 1), str(
                i)))
        os.makedirs(op.join(root, wfs_opt_path, str(iter_step + 1), str(i)))
        # prepare MOs in current step as MOs in next wfs opt step
        shutil.copy(
            op.join(root, wfs_opt_path, str(iter_step), str(i),
                    inp_file.split('.')[0] + ".gbw"),
            op.join(root, wfs_opt_path, str(iter_step + 1), str(i),
                    inp_file.split('.')[0] + "_last_step.gbw"))
        shutil.copy(
            op.join(root, wfs_opt_path, str(iter_step), str(i), inp_file),
            op.join(root, wfs_opt_path, str(iter_step + 1), str(i), inp_file))


def wfs_opt_further(iter_step):
    #generate M_values.dat, engrad.dat, energy.dat
    if (os.path.isdir(op.join(root, spin_ladder_path, str(iter_step)))):
        shutil.rmtree(op.join(root, spin_ladder_path, str(iter_step)))
    os.makedirs(op.join(root, spin_ladder_path, str(iter_step)))
    parse_parameters(iter_step)
    write_energies(iter_step)
    write_engrad(iter_step)
    write_S2(iter_step)
    write_spin_moment(iter_step)


def parse_parameters(iter_step):
    global orca_out_files
    global orca_engrad_files
    global N_Spin
    global N_BS
    global N_Engrad

    orca_out_files = []
    orca_engrad_files = []
    orca_engrad_files = sorted(orca_engrad_files_custom or engrad_files_list(
        op.join(root, wfs_opt_path, str(iter_step))))
    # clear the output files information from last step
    orca_out_files = sorted(orca_out_files_custom or output_files_list(op.join(
        root, wfs_opt_path, str(iter_step))))
    N_Engrad = len(parseEngrad(orca_engrad_files[0]))
    N_Spin = len(parseSpinMoment(orca_out_files[0]))
    N_BS = len(orca_out_files)


def orca_finished_succeed(fullname):
    orca_out = delegator.run('tail ' + '-n 2 ' + fullname + '| head -n 1 ').out #changed '.std_out' in '.out'
    if orca_finished_keywords in orca_out:
        return True
    else:
        return False


def geo_updated(fullname):
    if (time.time() - os.path.getmtime(fullname)) < waiting_time_wfs_finished:
        return True
    else:
        return False

def submit_orca_wfs(iter_step):
    node_name = open(op.join(root, 'machinefile')).readlines()
    for i in sorted(spin_states):
        run_dir =op.join(root, wfs_opt_path, str(iter_step), str(i)) 
        os.chdir(run_dir)
       #clean output files from last steps
        if inp_file not in os.listdir('.'):
            print("Oops! There is not input file for wfs opt")
        if clean_output and out_file in os.listdir('.'):
            os.remove(out_file)
        if (i=='0' or i=='6'): 
                orca_wfs_run_1 = delegator.run("ssh " + node_name[0].rstrip() + " 'cd " + run_dir + "; " + ORCA_EXE + " " + inp_file + " -d -v > " + out_file + " 2> " + err_file + "'", block=False)
        elif (i=='1' or i=='2'):
                orca_wfs_run_2 = delegator.run("ssh " + node_name[1].rstrip() + " 'cd " + run_dir + "; " + ORCA_EXE + " " + inp_file + " -d -v > " + out_file + " 2> " + err_file + "'", block=False)
        elif (i=='3' or i=='4'):
                orca_wfs_run_3 = delegator.run("ssh " + node_name[2].rstrip() + " 'cd " + run_dir + "; " + ORCA_EXE + " " + inp_file + " -d -v > " + out_file + " 2> " + err_file + "'", block=False)
        elif (i=='5' or i=='7'):
                orca_wfs_run_4 = delegator.run("ssh " + node_name[3].rstrip() + " 'cd " + run_dir + "; " + ORCA_EXE + " " + inp_file + " -d -v > " + out_file + " 2> " + err_file + "'", block=False)


# M_values.dat
def parseSpinMoment(orca_out_file):
    orca_output_stream = open(orca_out_file, "r").read()
    breakLine = Literal('\n')
    mullikenStartLine = "MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS"
    mullikenEndLine = "--------------------------------------------"
    point = Literal('.')
    e = CaselessLiteral('E')
    plusorminus = Literal('+') | Literal('-')
    number = Word(nums)
    integer = Combine(Optional(plusorminus) + number)
    floatNumber = Combine(integer + point + Optional(number) + Optional(
        e + integer))
    charge_spin = OneOrMore(floatNumber)
    Multinteger = OneOrMore(integer)
    pattern_mulliken = mullikenStartLine + LineEnd(
    ) + mullikenEndLine + LineEnd() + SkipTo(mullikenEndLine)
    mulliken_section = pattern_mulliken.searchString(
        orca_output_stream).asList()[-1][-1]
    pattern_spin_center = integer + SpinCenter + ":" + charge_spin + LineEnd()
    mulliken_spin_center = pattern_spin_center.searchString(
        mulliken_section).asList()
    mulliken_spin_population = []
    for i in mulliken_spin_center:
        mulliken_spin_population.append(i[4])
    return mulliken_spin_population


# S2_tot.dat
def parse_spin_square(orca_out_file):
    orca_output_stream = open(orca_out_file, "r").read()
    energyStartLine = "Expectation value of <S**2>     :"
    point = Literal('.')
    e = CaselessLiteral('E')
    plusorminus = Literal('+') | Literal('-')
    number = Word(nums)
    integer = Combine(Optional(plusorminus) + number)
    floatNumber = Combine(integer + point + Optional(number) + Optional(
        e + integer))
    pattern_en = energyStartLine + floatNumber + LineEnd()
    pattern_en.searchString(orca_output_stream).asList()[-1][1]
    return pattern_en.searchString(orca_output_stream).asList()[-1][1].split(
        '\n')

#engrad.dat


def parseEngrad(orca_engrad_file):
    engrad_output_stream = open(orca_engrad_file, "r").read()
    breakLine = Literal('\n')
    engradStartLine = "# The current gradient in Eh/bohr"
    keywordMeta = Literal('#')
    nParas = Word(nums)
    ntrunc = 3
    point = Literal('.')
    e = CaselessLiteral('E')
    plusorminus = Literal('+') | Literal('-')
    number = Word(nums)
    integer = Combine(Optional(plusorminus) + number)
    floatNumber = Combine(integer + point + Optional(number) + Optional(
        e + integer))
    engrad_paras = OneOrMore(floatNumber)
    Multinteger = OneOrMore(integer)
    pattern_engrad = engradStartLine + LineEnd() + keywordMeta + SkipTo(
        keywordMeta)
    return pattern_engrad.searchString(engrad_output_stream).asList()[0][
        -1].split('\n  ')


#Energies.dat
def parseEn(orca_out_file):
    orca_output = open(orca_out_file)
    orca_output_lines = orca_output.readlines()
    orca_output.close()
    for i in reversed(orca_output_lines):
        if "FINAL SINGLE POINT ENERGY" in i:
            return i.split()[-1]


def write_energies(iter_step):
    thefile = open(
        op.join(root, spin_ladder_path, str(iter_step), 'Energies.dat'), 'w')
    for orca_out_file in orca_out_files:
        thefile.write(parseEn(orca_out_file))
        thefile.write("\n")
    thefile.close()


def write_engrad(iter_step):
    thefile = open(
        op.join(root, spin_ladder_path, str(iter_step), 'engrad.dat'), 'w')
    thefile.write("%s " % N_BS)
    thefile.write("%s \n" % N_Engrad)
    for orca_engrad_file in orca_engrad_files:
        for i in parseEngrad(orca_engrad_file):
            thefile.write("%s " % i)
    thefile.close()


def write_S2(iter_step):
    thefile = open(
        op.join(root, spin_ladder_path, str(iter_step), 'S2_tot.dat'), 'w')
    for orca_out_file in orca_out_files:
        for i in parse_spin_square(orca_out_file):
            thefile.write("%s " % i)
            thefile.write("\n")
    thefile.close()


def write_spin_moment(iter_step):
    thefile = open(
        op.join(root, spin_ladder_path, str(iter_step), 'M_values.dat'), 'w')
    thefile.write("%s " % N_BS)
    thefile.write("%s \n" % N_Spin)
    for orca_out_file in orca_out_files:
        for i in parseSpinMoment(orca_out_file):
            thefile.write("%s " % i)
        thefile.write("\n")
    thefile.close()


def engrad_files_list(wfs_path):
    for dirpath, dirnames, filenames in os.walk(wfs_path):
        for filename in [f for f in filenames if f.endswith(".engrad")]:
            orca_engrad_files.append(os.path.join(dirpath, filename))
    return orca_engrad_files


def output_files_list(wfs_path):
    for dirpath, dirnames, filenames in os.walk(wfs_path):
        for filename in [f for f in filenames if f.endswith(".out")]:
            orca_out_files.append(os.path.join(dirpath, filename))
    return orca_out_files


if __name__ == "__main__":
    main()
