# copy the necessary files into this director

"""
Energy.cpp  Energy.o   Ewald.h  Field.cpp  Field.o         NbrListPBC.h  Potential.h  TwoBody.h
Energy.h    Ewald.cpp  Ewald.o  Field.h    NbrListPBC.cpp  NbrListPBC.o  TwoBody.cpp  TwoBody.o
"""
import os
import sys
import subprocess
import shutil

files = []
files.append("Energy.cpp")
files.append("Energy.h")
files.append("Ewald.cpp")
files.append("Ewald.h")
files.append("Field.cpp")
files.append("Field.h")
files.append("Potential.h")
files.append("TwoBody.cpp")
files.append("TwoBody.h")
files.append("NbrListPBC.cpp")
files.append("NbrListPBC.h")

path_dir = "/mnt/Progs/workspace/rdforce/source"

new_dir = "pair"

if os.path.exists(new_dir) != True:
    os.mkdir(new_dir)

for f in files:
    source = os.path.join(path_dir, f)
    dest = os.path.join(new_dir, f)

    cmd = "cp " + source + " " + dest

    handle = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    stdout, stderr = handle.communicate()

    # Format stdout to utf-8 for python 3, python 2 should be untouched.
    if not isinstance(stdout, str):
        stdout = stdout.decode("utf-8")

    # Format stderr to utf-8 for python 3, python 2 should be untouched.
    if not isinstance(stderr, str):
        stderr = stderr.decode("utf-8")

    # Grab the return code.
    errorstate = handle.returncode

    #submit first line
    if errorstate != 0:
        print(stderr)
        sys.exit()

print("finished copying pair-potential files ")
