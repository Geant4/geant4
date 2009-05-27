#!/usr/bin/env python

import os, sys
import datetime

#vfem = os.environ.get("VFEM")
#dr = vfem + "/test44"
#os.chdir(dr)

#REFERENCE = datetime.datetime.now().strftime('%d_%m_%Y-%H:%M:%S')
REFERENCE = datetime.datetime.now().strftime('%d_%m_%Y')

if not os.path.isdir(REFERENCE):
        os.mkdir(REFERENCE)

os.chdir(REFERENCE)

for file in os.listdir(os.getcwd()):
    if (os.path.splitext(file)[1] == ".txt"):
        os.remove(file)

g4ins = os.environ.get("G4INSTALL")
dir  = g4ins + "/tests/test44/"
src  = dir + "Exp_Data/"
dst = os.getcwd()
names = os.listdir(src)
for name in names:
    if (os.path.splitext(name)[1] == ".txt"):
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        os.symlink(srcname, dstname)

PHYSLIST = "QBBC"

wat = "_water_"
tPart = "p", "he4", "c12"
phys = "opt0", "opt2", "opt3"

g4bin = os.environ.get("G4BIN")
g4sys = os.environ.get("G4SYSTEM")
work = g4bin + "/" + g4sys + "/test44"

for ph in phys:
    f1 = wat + ph
    f2 = f1 + ".log"
    for pp in tPart:
        f = pp + f2
        try:
            # remove old file, if any
            os.remove(f)
        except os.error:
            pass
        cmd = work + " " + dir + pp + f1 + ".in >& " + f
        os.system(cmd)
        snam = "Bragg.out"
        dnam = pp + "_" + ph + ".out"
        os.rename(snam, dnam)

pyfile = dir + "utils/reader_test44-online.py"
#os.chmod(pyfile, 755)

global tyPart, refer
refer = sys.argv[1]
for pp in tPart:
    tyPart = pp
    execfile(pyfile, globals())


