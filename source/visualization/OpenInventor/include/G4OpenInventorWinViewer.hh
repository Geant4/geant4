//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4OpenInventorWinViewer.hh 66373 2012-12-18 09:41:34Z gcosmo $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.

#ifndef G4OPENINVENTORWINVIEWER_HH
#define G4OPENINVENTORWINVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4OpenInventorViewer.hh"

#include <windows.h>

class Geant4_SoWinExaminerViewer;

class G4OpenInventorWinViewer: public G4OpenInventorViewer {
public: //G4VViewer
  virtual void FinishView();
  virtual void SetView();
protected:
  virtual void ViewerRender();
  virtual SoCamera* GetCamera();
public:
  G4OpenInventorWinViewer(G4OpenInventorSceneHandler& scene,
		       const G4String& name = "");
  virtual ~G4OpenInventorWinViewer();
  void Initialise();

private:
  static LRESULT CALLBACK WindowProc(HWND,UINT,WPARAM,LPARAM);
private:
  HWND fShell;
  Geant4_SoWinExaminerViewer* fViewer;
};

#endif

#endif
