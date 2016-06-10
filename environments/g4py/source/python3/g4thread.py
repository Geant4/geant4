"""
# ==================================================================
#   Python module
#
#   Geant4 threading module
#
#                                              Q, 2005
# ==================================================================
"""
#$Id: g4thread.py 66892 2013-01-17 10:57:59Z gunter $

import thread
from G4run import *

# ------------------------------------------------------------------
# BeamOn in a new thread
# ------------------------------------------------------------------
def _TBeamOn(self, nevent):
  "generate events in a thread"
  args = (nevent,)
  thread.start_new_thread(self.BeamOn, args)

G4RunManager.TBeamOn= _TBeamOn

