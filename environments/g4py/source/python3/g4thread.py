"""
# ==================================================================
#   Python module
#
#   Geant4 threading module
#
#                                              Q, 2005
# ==================================================================
"""

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

