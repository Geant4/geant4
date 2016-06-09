#$Id: g4viscp.py,v 1.1 2010-12-02 08:22:21 kmura Exp $
"""
# ==================================================================
#   Python module
#
#   Visualization Control Panel
#
#                                              Q, 2005
# ==================================================================
"""
from G4interface import *

# ------------------------------------------------------------------
# Scene
# ------------------------------------------------------------------
class G4Scene :
  "Scene"
  def __init__(self, aname, vol= "world", acopyno=0,
               amode=0, bmode=1):
    self.name=    aname
    self.volume=  vol
    self.copyno=  acopyno
    self.mode_eventaction= amode  # 0: accumulate / 1: refresh
    self.mode_runaction=   bmode  # 0: accumulate / 1: refresh
    self.mode=    ("accumulate", "refresh")

  def create_scene(self):
    ApplyUICommand("/vis/scene/create " + self.name)
    ApplyUICommand("/vis/scene/add/volume %s %d" %
                   (self.volume, self.copyno))
    ApplyUICommand("/vis/scene/add/trajectories")
    self.update_scene()
      
  def update_scene(self):
    ApplyUICommand("/vis/scene/select " + self.name)
    ApplyUICommand("/vis/sceneHandler/attach")
    ApplyUICommand("/vis/scene/endOfEventAction %s" %
                   (self.mode[self.mode_eventaction]) )
    ApplyUICommand("/vis/scene/endOfRunAction %s" %
                   (self.mode[self.mode_runaction]) )

# ------------------------------------------------------------------
# Visualization Control Panel
# ------------------------------------------------------------------
class G4VisCP :
  "G4 Visualization Control Panel"

  def __init__(self, gsys="OGLIX"):
    self.gsystem=    gsys    
    self.scenelist=  [G4Scene("default")]    
    self.viewpoint=  [270., 90.]

    rc= ApplyUICommand("/vis/open " + gsys)
    if (rc != 0):
      return

    self.scenelist[0].create_scene()
    ApplyUICommand("/vis/viewer/set/viewpointThetaPhi %f %f"
                   % (self.viewpoint[0], self.viewpoint[1]) )
    ApplyUICommand("/tracking/storeTrajectory 1")

  def add_scene(self, ascene):
    self.scenelist.append(ascene)
    
  def select_scene(self, iscene):
    self.scenelist[iscene].update_scene()
    ApplyUICommand("/vis/viewer/set/viewpointThetaPhi %f %f"
                   % (self.viewpoint[0], self.viewpoint[1]) )

