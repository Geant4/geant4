"""
# ==================================================================
#   Python module
#
#   Geant4 threading module
#
#                                              Q, 2005
# ==================================================================
"""
#$Id: g4thread.py,v 1.3 2008-12-03 07:01:04 kmura Exp $

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

