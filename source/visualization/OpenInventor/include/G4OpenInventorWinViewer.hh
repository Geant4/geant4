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
// $Id: G4OpenInventorWinViewer.hh,v 1.7 2004/11/25 13:39:54 gbarrand Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
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
protected:
  virtual void ViewerRender();
  virtual SoCamera* GetCamera();
public:
  G4OpenInventorWinViewer(G4OpenInventorSceneHandler& scene,
		       const G4String& name = "");
  virtual ~G4OpenInventorWinViewer();
private:
  static LRESULT CALLBACK WindowProc(HWND,UINT,WPARAM,LPARAM);
private:
  HWND fShell;
  Geant4_SoWinExaminerViewer* fViewer;
};

#endif

#endif
