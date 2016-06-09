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
// $Id: G4OpenInventorWinViewer.hh,v 1.2 2004/04/08 10:49:10 gbarrand Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORWINVIEWER_HH
#define G4OPENINVENTORWINVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include "G4VViewer.hh"

#include <windows.h>

class SoWinExaminerViewer;
class SoSelection;
class G4OpenInventorSceneHandler;
class G4VInteractorManager;

class G4OpenInventorWinViewer: public G4VViewer {

public:
  G4OpenInventorWinViewer(G4OpenInventorSceneHandler& scene,
		       const G4String& name = "");
  virtual ~G4OpenInventorWinViewer();
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
  HWND fShell;
  SoWinExaminerViewer* fViewer;
  SoSelection* fSelection;
  G4VInteractorManager* fInteractorManager;
};

#endif

#endif
