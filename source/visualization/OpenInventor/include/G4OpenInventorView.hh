// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorView.hh,v 1.1 1999-01-07 16:15:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor view - opens window, hard copy, etc.

#ifndef G4OPENINVENTORVIEW_HH
#define G4OPENINVENTORVIEW_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include <rw/tvordvec.h>

#include <Inventor/Xt/SoXt.h>
#include <Inventor/Xt/SoXtComponent.h>
#include <Inventor/nodes/SoSelection.h>

#include "G4VView.hh"

class G4OpenInventorScene;
class G4VInteractorManager;

//
// Base class for various OpenInventorView classes.
//
class G4OpenInventorView: public G4VView {

public:
  		  G4OpenInventorView (G4OpenInventorScene& scene,
				      const G4String& name = "");
  		 ~G4OpenInventorView ();
  void 		  DrawView ();
  void 		  ShowView ();
  G4bool 	  GetOIVisualFound () const;

private:
  void   	  ClearView           ();
  void   	  FinishView          ();
  void   	  SetView             ();
  void   	  KernelVisitDecision ();

  G4OpenInventorScene&	fScene; 	// Graphics Scene for this view.
  G4bool		OIvisualfound;
 
  Widget 	        G4OIShell;
  SoXtComponent*        G4OIViewer;	  // The Inventor Viewer
  SoSelection*          G4OISelection;
  G4VInteractorManager* interactorManager;
};

inline G4bool G4OpenInventorView::GetOIVisualFound () const {
  return OIvisualfound;
}

#endif

#endif


