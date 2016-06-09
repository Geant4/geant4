"""
# ==================================================================
#   Python module
#
#   Geant4 threading module
#
#                                              Q, 2005
# ==================================================================
"""
#$Id: g4thread.py,v 1.1 2010-12-02 08:22:21 kmura Exp $

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

