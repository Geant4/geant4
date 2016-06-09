# $Id: liz.py,v 1.2 2003/06/16 17:06:44 dressel Exp $
# -------------------------------------------------------------------
# GEANT4 tag $Name: geant4-06-00 $
# -------------------------------------------------------------------
#

import os, sys, string

if ( not os.environ.has_key("ANAPHETOP") ) :
    os.environ["ANAPHETOP"] = "/afs/cern.ch/sw/lhcxx"


if ( not os.environ.has_key("PLATF") ) :
    os.environ["PLATF"] = "redhat73/gcc-3.2"

if ( not os.environ.has_key("ANAPHE_VERSION") ) :
    if ( not os.environ.has_key("ANAPHEVERS") ) :
        os.environ["ANAPHEVERS"] = "5.0.4"
else :
    os.environ["ANAPHEVERS"] = os.environ["ANAPHE_VERSION"]
    
if ( not os.environ.has_key("PUBDOMVERS") ) :
    os.environ["PUBDOMVERS"] = "2.0.0"

if ( not os.environ.has_key("ANAPHE_REL_DIR") ) :
    os.environ["ANAPHE_REL_DIR"] = os.environ["ANAPHETOP"]+"/specific/"+os.environ["PLATF"]+"/"+os.environ["ANAPHEVERS"]


if ( not os.environ.has_key("OS") ) :
    os.environ["OS"] = "Linux"  # for now !!!

# for debugging purposes this might be set differently ...
if ( not os.environ.has_key("LIZARD_ROOT") ) :
    os.environ["LIZARD_ROOT"] = os.environ["ANAPHE_REL_DIR"] + "/python"

os.environ["LIZARD_LIB"] = os.environ["LIZARD_ROOT"] + "/lib"

#-toDo: clean up the LD_LIBRARY_PATH, PATH and PYTHONPATH variables from wrong versions 
# of Anaphe s/w
#export LD_LIBRARY_PATH=`${LIZARD_ROOT}/bin/cleanupPath.py LD_LIBRARY_PATH`
#export PATH=`${LIZARD_ROOT}/bin/cleanupPath.py PATH`
#export PYTHONPATH=`${LIZARD_ROOT}/bin/cleanupPath.py PYTHONPATH`

# add to path in order to find xmgrace ...
if (not os.environ.has_key("GRACE_DIR") ) :
    os.environ["GRACE_DIR"] = os.environ["ANAPHETOP"]+"/specific/"+os.environ["PLATF"]+"/PublicDomainPackages/" + os.environ["PUBDOMVERS"] + "/grace"

os.environ["PATH"] = os.environ["GRACE_DIR"] + "/bin:" + os.environ["PATH"]
os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":" + os.environ["LIZARD_LIB"] + ":" + os.environ["ANAPHE_REL_DIR"] + "/lib"


# toDo:                                
# if [ `uname` = "SunOS" ] ; then
#   SUN_CC_DIR=/afs/cern.ch/project/sun/solaris/opt/SUNWspro62Apr02
#   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${SUN_CC_DIR}/lib:/usr/local/lib   # /opt/SUNWspro/lib:/usr/local/lib:/opt/SUNWspro/WS6U1/lib
# fi
# 
# to find the proper versions of libs for python-2.0 and swig-1.3a5 
# (it should be SWIG_DIR)

if (not os.environ.has_key("SWIG_DIR")) :
    os.environ["SWIG_DIR"] = os.environ["ANAPHETOP"]+"/specific/"+os.environ["PLATF"]+"/PublicDomainPackages/" + os.environ["PUBDOMVERS"]


os.environ["LD_LIBRARY_PATH"] = os.environ["SWIG_DIR"] + "/lib:" + os.environ["LD_LIBRARY_PATH"]
    
# need to add "." for automatic code compilation (ntuples et al)
os.environ["LD_LIBRARY_PATH"] = os.environ["LD_LIBRARY_PATH"] + ":" + os.environ["LIZARD_ROOT"] + "/" + os.environ["PLATF"] + ":."

# add to PYTHONPATH
sys.path.append(os.environ["LIZARD_LIB"])
sys.path.append(os.environ["LIZARD_ROOT"] + "/" + os.environ["PLATF"])
sys.path.append(os.environ["LIZARD_ROOT"] + "/src")
sys.path.append(os.environ["HOME"])

sys.path.append(os.environ["LIZARD_ROOT"] + "/contrib")

# toDo
#if [ -z $LIZARD_NO_PUBLIC_CONTRIB ] ; then
#  export PYTHONPATH=${PYTHONPATH}:${ANAPHETOP}/share/PythonContrib
#fi

# # define location of python
# if [ -z $PYDIR ] ; then
#   export PYDIR=${ANAPHETOP}/specific/${PLATF}/PublicDomainPackages/${PUBDOMVERS}
# fi
# 
# if [ -z $PYTHON ] ; then 
#   export PYTHON=${PYDIR}/bin/python2.2
# fi
# 

# create a local temporary file (with the PID as part of the name) to 
# speed up access from Nag fitter if starting directory is in AFS:
if (not os.environ.has_key("LIZARD_KEEP_FITRESULT")) :
   os.system("rm -f e04ucc.r")
   os.system("ln -s /tmp/pid-$$-e04ucc.r e04ucc.r")

# $PYTHON -i ${LIZARD_ROOT}/bin/.Lizardrc $* 
print "executing startup from'"+os.environ["LIZARD_ROOT"] + "/bin/.Lizardrc"+"'"
# execfile(os.environ["LIZARD_ROOT"] + "/bin/.Lizardrc")



# Settings for Lizard 

lizardVersion = "3.0.0.5"
lizardVersionDate = "20 Dec 2002"

lizardBatch  = None
lizardObjy   = None
lizardNag    = None
lizardNoGraphics = None

from math import *
from time import *

import string
import os
import sys

import random
import atexit
import gc 

# the following works only for python-2.0 (or later)

# for command-line completion :
import rlcompleter
rlcompleter.readline.parse_and_bind("tab: complete")

try :
   if (sys.version_info[0] == 2) :
      # for history (across-sessions)
      import readline
      # limit history file to 1000 lines
      readline.set_history_length(1000) 
      # read previous file (if existing)
      histfile = os.path.join(os.environ["HOME"], ".LizHist")
      try:
        readline.read_history_file(histfile)
      except IOError:
        print "cannot read command history file ", histfile
        pass
      # register function to write history file at exit
      atexit.register(readline.write_history_file, histfile)
      del histfile
      # end history
except :
   print "\nerror while accessing command history file"
   raise

# end python-2.0 specific part ...


# --------------------------------------------------------------------------------

def usage ():
   print """

usage: startLizard.sh [<options>] [<file1> <file2> ...]

where the optional <options> can be one of the following:

       -?, -h, --help : print this text
       -v, --version  : print version info

       -b, --batch    : run in batch mode and execute the scripts passed after all options

       --noGraphics   : don't instantiate a Plotter (use this if your DISPLAY variable is not set)

       --useNag       : use minimizer engine from Nag-C library

       --useObjy      : use Objectivity/DB for persistent Histograms and (row-wise) Ntuples

         default      : use FML with Minuit minimizer engine

the (optional) list of files will be executed after startup

NOTE: the options need to be _before_ any script file

"""


# --------------------------------------------------------------------------------
# set some defaults: use FML, Objy, Nag-C
# --------------------------------------------------------------------------------
lizardFML    = 1
lizardHBook  = 0
lizardNag    = 0
# --------------------------------------------------------------------------------

# --------------------------------------------------------------------------------

def lizardVersionInfo() :
   print "\nThis is Lizard version " + lizardVersion + "\n"


# --------------------------------------------------------------------------------
# check if we got any flags:
# --------------------------------------------------------------------------------

import getopt
#optlist = []
#args = []
try:
   optlist, args = getopt.getopt(sys.argv[1:], ['?', 'h', 'v', 'b'], 
                                 ['help', 'version', 'batch', 'noGraphics',
                                  'useNag', 'useNag', 'useObjy', 'objy' ])

except :
   print "\nunknown option:",o,"\n"
   usage()
   sys.exit()    

# --------------------------------------------------------------------------------
# have a first look at options now, see if the user wants some info without
# the need to start the system ...
# --------------------------------------------------------------------------------

for o, a in optlist:
    if o in ("-?", "-h", "--help"):
       usage()
       sys.exit()
    elif o in ("-v", "--version",):
       lizardVersionInfo()
       sys.exit()
    elif o in ("", "--useObjy", "--objy",):
       lizardObjy = 1


import dl 
# --------------------------------------------------------------------------------
# check if Objectivity can be used, if requested
# --------------------------------------------------------------------------------

if (lizardObjy == 1) :
  try:
    dl.open( "liboo.so" )
  except:
    print "Objectivity requested, but liboo not found in LD_LIBRARY_PATH !! Aborting !"
    sys.exit()

# --------------------------------------------------------------------------------
# check if Nag_C can be used, 
# --------------------------------------------------------------------------------

if (lizardNag == 1 ) :
   try:
      dl.open( "libnagc.so" )
   except:
      print "Minimizer from NAG requested, but libnagc not found in LD_LIBRARY_PATH !! Aborting !"
      sys.exit()

# --------------------------------------------------------------------------------
# check the other possible flags:
# --------------------------------------------------------------------------------
for o, a in optlist:
    if o in ("-b", "--batch",):
       lizardBatch = 1
    elif o in ("--noGraphics",):
       lizardNoGraphics = 1
    elif o in ("--nag", "--useNag"):
       lizardNag = 1

# --------------------------------------------------------------------------------

print "\n"

# --------------------------------------------------------------------------------
# check now whether DISPLAY is set.
# --------------------------------------------------------------------------------
if (lizardNoGraphics != 1) :
   try:
     disp = os.environ["DISPLAY"]
   except:
     print "\n==========> no DISPLAY set, switching to no-graphics mode."
     lizardNoGraphics = 1

# prompt
sys.ps1=":-) "

# --------------------------------------------------------------------------------

# welcome message parameters
lizMsgWelcome     = "Welcome to Lizard"
lizMsgID          = "Version "+ lizardVersion + " (" + lizardVersionDate + ")"
lizMsgURL         = "http://cern.ch/Anaphe"
lizMsgSide        = "|"
lizMsgCorner      = "+"
lizMsgTop         = "-"
lizMsgWidth       = 50
lizMsgSideIndent  =  8
lizMsgTopIndent   =  1

# generate welcome message
lizMsgBar     = ""
for i in range (0,lizMsgWidth) :
    lizMsgBar = lizMsgBar + lizMsgTop

lizMsgPad     = ""
for i in range (0,lizMsgSideIndent) :
    lizMsgPad = lizMsgPad + " "

lizMsgBlank = lizMsgSide + string.center("",lizMsgWidth) + lizMsgSide
for i in range (0, lizMsgTopIndent) :
    print ""
    
print lizMsgPad + lizMsgCorner + lizMsgBar + lizMsgCorner
print lizMsgPad + lizMsgBlank 
print lizMsgPad + lizMsgSide + string.center(lizMsgWelcome,lizMsgWidth) + lizMsgSide
print lizMsgPad + lizMsgBlank
print lizMsgPad + lizMsgSide + string.center(lizMsgID,lizMsgWidth) + lizMsgSide
print lizMsgPad + lizMsgBlank
print lizMsgPad + lizMsgSide + string.center(lizMsgURL,lizMsgWidth) + lizMsgSide
print lizMsgPad + lizMsgBlank
print lizMsgPad + lizMsgCorner + lizMsgBar + lizMsgCorner
for i in range (0, lizMsgTopIndent) :
    print ""


# --------------------------------------------------------------------------------

print "Chosen configuration: ",
if ( lizardObjy == 1 ) : 
   print "Objectivity",
else:
   print "HBook and XML",
print " persistency and ",
if ( lizardNag == 1 ) : 
   print "Nag-C",
else:
   print "Minuit",
print "minimizer engine"

print "\nType help() for help\n"

sys.stdout.flush()

# --------------------------------------------------------------------------------
# global startup files ...
# --------------------------------------------------------------------------------

# remember where we started from ...
curDir = os.getcwd()

# check for global StartupFiles (so far, the order doesn't matter):
startDir = os.environ['LIZARD_ROOT']+"/startUpFiles/"
fileList = os.listdir(startDir)

# change to directory containing the startup files and exec them one by one:
os.chdir(startDir)
try:
  for file in fileList :
    if (string.find(file,".py") == -1) :
       continue;
    if (file == "constants.py") :
       continue;
    if (os.path.isfile(file)) :
        # print "loading startUp file ", file
        execfile(file)
    else :
        print "Strange: cannot \"execfile\" startUp file ", file, " in ", startDir, " !?!?"
except:
   raise

#-ap execfile("constants.py")
# return to where we left from ...
os.chdir(curDir)

LizardIsInitialized = 1

print "\n\nLizard initialised\n"

# --------------------------------------------------------------------------------
# local/user-specific startup files ...
# --------------------------------------------------------------------------------

# check for local initfiles:
homeStart = os.environ['HOME']+"/.Lizardrc"
if os.path.isfile(homeStart):
   execfile(homeStart)

# avoid reading it twice ...
if (os.getcwd() != os.environ['HOME'] and 
    os.getcwd() != os.environ['LIZARD_ROOT']+"/bin" ) :
    localStart = "./.Lizardrc"	
    if os.path.isfile(localStart):
       execfile(localStart)

# --------------------------------------------------------------------------------
def lizardCleanUp() :
  global af, tf, pf, pl
  # in reverse order of construction 
  if (lizardNoGraphics != 1) :
    pl.thisown = 1 ; del pl
    pf.thisown = 1 ; del pf
  tf.thisown = 1 ; del tf
  af.thisown = 1 ; del af
  # print "global objects deleted ... " ; sys.stdout.flush()

  for lib in theListOfLoadedLibraries.keys() :
     sys.stdout.flush()
     theListOfLoadedLibraries[lib].close()
  # print "\nall libs unloaded ... " ; sys.stdout.flush()

  print "\nobjects deleted and libs unloaded, exiting python ... " ; sys.stdout.flush()

  return
# --------------------------------------------------------------------------------
def exit() :
  lizardCleanUp()
  sys.exit()

# --------------------------------------------------------------------------------
# execute file given at input (ignoring startup file .Lizardrc) ...
# --------------------------------------------------------------------------------

if (len(args) > 0 ) :
  args[0] = os.path.abspath(args[0])
  file = args[0]
  # reset the environment ... the scripts will now look as if invoked from the shell ...
  sys.argv = args
  print "executing ", file
  if (os.path.isfile(file)) :
    try :
      execfile(file)
    except:
      # if something goes wrong in batch mode, clean up and exit ... 
      if (lizardBatch == 1) : 
  	lizardCleanUp()
  	sys.exit(2)
      else :
  	raise
  else :
      print "file ", file, " not found."


if (lizardBatch == 1) :
    # clean up (explicitly destroy some objects)
    lizardCleanUp()
    sys.exit()



