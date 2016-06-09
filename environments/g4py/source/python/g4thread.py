#$Id: g4thread.py,v 1.2 2006/04/25 08:09:45 kmura Exp $
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
def TBeamOn(self, nevent):
  "generate events in a thread"
  args= (nevent,)
  thread.start_new_thread(self.BeamOn, args)

G4RunManager.TBeamOn= TBeamOn

