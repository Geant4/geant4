// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorViewer.hh,v 1.4 1999-05-12 14:00:45 barrand Exp $
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

#include "G4VViewer.hh"

class SoXtComponent;
class SoSelection;
class SoWindow;
class G4OpenInventorSceneHandler;
class G4VInteractorManager;

//
// Base class for various OpenInventorView classes.
//
class G4OpenInventorViewer: public G4VViewer {

public:
  G4OpenInventorViewer (G4OpenInventorSceneHandler& scene,
			const G4String& name = "");
  virtual ~G4OpenInventorViewer ();
  void DrawView ();
  void ShowView ();
private:
  void ClearView           ();
  void FinishView          ();
  void SetView             ();
  void KernelVisitDecision ();

  G4OpenInventorSceneHandler&	fSceneHandler; 	// Graphics Scene for this view.
  Widget fShell;
  SoWindow* fWindow;
  SoXtComponent* fViewer;	  // The Inventor Viewer
  SoSelection* fSelection;
  G4VInteractorManager* fInteractorManager;
};

#endif

#endif


