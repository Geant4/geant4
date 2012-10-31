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
// $Id: G4OpenInventorXtViewer.hh,v 1.15 2009-09-18 12:48:43 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Jeff Kallenbach 01 Aug 1996
// OpenInventor viewer - opens window, hard copy, etc.
// Frederick Jones and TJR  October 2012
// Extended viewer based on G4OpenInventorXtViewer.hh
// Uses G4OpenInventorXtExaminerViewer.

#ifndef G4OPENINVENTORXTEXTENDEDVIEWER_HH
#define G4OPENINVENTORXTEXTENDEDVIEWER_HH

#ifdef G4VIS_BUILD_OI_DRIVER

// Inheritance :
#include "G4OpenInventorXtViewer.hh"

#include <X11/Intrinsic.h>
#include <Inventor/nodes/SoEventCallback.h>

class G4OpenInventorXtExaminerViewer;

class G4OpenInventorXtExtendedViewer: public G4OpenInventorXtViewer {
public:
  G4OpenInventorXtExtendedViewer(G4OpenInventorSceneHandler& scene,
		         const G4String& name = "");
  void Initialise();
  virtual ~G4OpenInventorXtExtendedViewer();
protected:
  static void EscapeFromKeyboardCbk(void * o);
  G4OpenInventorXtExaminerViewer* fViewer;
};

#endif

#endif
