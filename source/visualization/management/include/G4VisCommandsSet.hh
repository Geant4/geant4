// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsSet.hh,v 1.6 2001-02-05 02:34:01 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// /vis~/set/ commands
// John Allison  7th September 1997

#ifndef G4VISCOMMANDSSET_HH
#define G4VISCOMMANDSSET_HH

#include "globals.hh"
#include "G4VisManager.hh"
#include "G4ViewParameters.hh"

////////////////////////////////////////////////////  /vis~/set/...  ////
//vis \hline
//vis /vis~/set/ &&
//vis ...menu of set commands. \\%
class G4VisCommandSet {
public:
  G4String GetCommandName () const {return "/vis~/set/";}
  G4String GetGuidance () const {
    return "...menu of set commands.";
  }
};

///////////////////////////////////////////  /vis~/set/culling  ////
//set \hline
//set /vis~/set/culling & true/false &
//set Global culling flag.  Does not change specific culling flags. \\%
class G4VisCommandSetCulling {
public:
  G4String GetCommandName () const {return "/vis~/set/culling";}
  G4String GetGuidance () const {
    return "Global culling flag.  Does not change specific culling flags.";
  }
  G4String GetValueName () const {return "global culling flag";}
  G4bool GetValue () const {
    return G4VisManager::GetInstance () -> GetCurrentViewParameters ().
      IsCulling ();
  }
  void SetValue (G4bool value) {
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/set/culling\".");
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    pVMan -> SetCurrentViewParameters ().SetCulling (value);
    G4VViewer* pView = pVMan -> GetCurrentViewer ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (pVMan -> GetCurrentViewParameters ());
    }
    G4cout << "Issue Draw or refresh to see effect." << G4endl;
  }
};

//////////////////////////////////////////  /vis~/set/cull_covered_daughters
//set \hline
//set /vis~/set/ cull\_covered\_daughters & true/false &
//set Cull (i.e., do not Draw) daughters covered by opaque mothers. \\%
class G4VisCommandSetCullCoveredDaughters {
public:
  G4String GetCommandName () const {return "/vis~/set/cull_covered_daughters";}
  G4String GetGuidance () const {
    return "Cull (i.e., do not Draw) daughters covered by opaque mothers.";
  }
  G4String GetValueName () const {return "culling covered daughters flag";}
  G4bool GetValue () const {
    return G4VisManager::GetInstance () -> GetCurrentViewParameters ().
      IsCullingCovered ();
  }
  void SetValue (G4bool value) {
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/set/culling\".");
    G4cout << "\nNote: this is only effective in surface drawing style,"
      "\nand then only if the volumes are visible and opaque, and then"
      "\nonly if no sections or cutways are in operation."
	 << G4endl;
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    pVMan -> SetCurrentViewParameters ().SetCullingCovered (value);
    G4VViewer* pView = pVMan -> GetCurrentViewer ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (pVMan -> GetCurrentViewParameters ());
    }
    G4cout << "Issue Draw or refresh to see effect." << G4endl;
  }
};

//////////////////////////////////////////  /vis~/set/cull_invisible_objects
//set \hline
//set /vis~/set/ cull\_invisible\_objects & true/false &
//set Cull (i.e., do not Draw) ``invisible'' objects. \\%
class G4VisCommandSetCullInvisible {
public:
  G4String GetCommandName () const {return "/vis~/set/cull_invisible_objects";}
  G4String GetGuidance () const {
    return "Cull (i.e., do not Draw) \"invisible\" objects.";
  }
  G4String GetValueName () const {return "culling invisible objects flag";}
  G4bool GetValue () const {
    return G4VisManager::GetInstance () -> GetCurrentViewParameters ().
      IsCullingInvisible ();
  }
  void SetValue (G4bool value) {
    G4VisManager::PrintCommandDeprecation("Use \"/vis/viewer/set/culling\".");
    G4VisManager* pVMan = G4VisManager::GetInstance ();
    pVMan -> SetCurrentViewParameters ().SetCullingInvisible (value);
    G4VViewer* pView = pVMan -> GetCurrentViewer ();
    if (pView) {
      // Copy current view parameters into current view.
      pView -> SetViewParameters (pVMan -> GetCurrentViewParameters ());
    }
    G4cout << "Issue Draw or refresh to see effect." << G4endl;
  }
};

#endif
