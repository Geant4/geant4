import sys
origArgv = sys.argv
sys.argv = ["","--noGraphics"]
from liz import *
sys.argv = origArgv
