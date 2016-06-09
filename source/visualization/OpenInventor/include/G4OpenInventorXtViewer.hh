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
// $Id: G4OpenInventorXtViewer.hh,v 1.1 2004/04/08 09:37:30 gbarrand Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORXTVIEWER_HH
#define G4OPENINVENTORXTVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include "G4VViewer.hh"

#include <X11/Intrinsic.h>
class SoXtExaminerViewer;

class SoSelection;
class G4OpenInventorSceneHandler;
class G4VInteractorManager;

//
// Base class for various OpenInventorView classes.
//
class G4OpenInventorXtViewer: public G4VViewer {

public:
  G4OpenInventorXtViewer(G4OpenInventorSceneHandler& scene,
		       const G4String& name = "");
  virtual ~G4OpenInventorXtViewer();
  void DrawView();
  void ShowView();
private:
  void ClearView();
  void FinishView();
  void SetView();
  void KernelVisitDecision();
  G4bool CompareForKernelVisit(G4ViewParameters&);
  G4OpenInventorSceneHandler& fG4OpenInventorSceneHandler;
  G4ViewParameters fLastVP;  // Memory for making kernel visit decisions.
  Widget fShell;
  SoXtExaminerViewer* fViewer;
  SoSelection* fSelection;
  G4VInteractorManager* fInteractorManager;
};

#endif

#endif
