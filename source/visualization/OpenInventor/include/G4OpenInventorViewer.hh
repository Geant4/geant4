//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4OpenInventorViewer.hh,v 1.9 2001-08-14 18:38:34 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORVIEWER_HH
#define G4OPENINVENTORVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

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
  G4bool CompareForKernelVisit(G4ViewParameters&);

  G4OpenInventorSceneHandler& fG4OpenInventorSceneHandler;
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
  Widget fShell;
  SoWindow* fWindow;
  SoXtComponent* fViewer;	  // The Inventor Viewer
  SoSelection* fSelection;
  G4VInteractorManager* fInteractorManager;
};

#endif

#endif
