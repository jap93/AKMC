# copy the necessary files into this director

"""
Energy.cpp  Energy.o   Field.h  ManyBody.cpp  MetalPotential.cpp  NbrListPBC.cpp  NbrListPBC.o
Energy.h    Field.cpp  Field.o  ManyBody.h    MetalPotential.h    NbrListPBC.h
"""
import os
import sys
import subprocess
import shutil

files = []
files.append("Energy.cpp")
files.append("Energy.h")
files.append("Field.cpp")
files.append("Field.h")
files.append("ManyBody.cpp")
files.append("ManyBody.h")
files.append("MetalPotential.cpp")
files.append("MetalPotential.h")
files.append("NbrListPBC.cpp")
files.append("NbrListPBC.h")

path_dir = "/mnt/Progs/workspace/rdforce/metal"

new_dir = "metal"

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

print("finished copying metal files ")
