// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorViewer.hh,v 1.2 1999-01-11 00:47:52 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORVIEWER_HH
#define G4OPENINVENTORVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include <rw/tvordvec.h>

#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/SoXtComponent.h>
#include <Inventor/nodes/SoSelection.h>

#include "G4VViewer.hh"

class G4OpenInventorSceneHandler;
class G4VInteractorManager;

//
// Base class for various OpenInventorView classes.
//
class G4OpenInventorViewer: public G4VViewer {

public:
  		  G4OpenInventorViewer (G4OpenInventorSceneHandler& scene,
				      const G4String& name = "");
  		 ~G4OpenInventorViewer ();
  void 		  DrawView ();
  void 		  ShowView ();
  G4bool 	  GetOIVisualFound () const;

private:
  void   	  ClearView           ();
  void   	  FinishView          ();
  void   	  SetView             ();
  void   	  KernelVisitDecision ();

  G4OpenInventorSceneHandler&	fSceneHandler; 	// Graphics Scene for this view.
  G4bool		OIvisualfound;
 
  Widget 	        G4OIShell;
  SoXtComponent*        G4OIViewer;	  // The Inventor Viewer
  SoSelection*          G4OISelection;
  G4VInteractorManager* interactorManager;
};

inline G4bool G4OpenInventorViewer::GetOIVisualFound () const {
  return OIvisualfound;
}

#endif

#endif


