#$Id: g4thread.py,v 1.1 2006-02-27 09:56:05 kmura Exp $
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

