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
// Satoshi Tanaka  31th May 2001
// A dummy viewer for GAGTreeSceneHandler.

#include "G4GAGTreeViewer.hh"

#include "G4ios.hh"
#include "g4std/strstream"

G4GAGTreeViewer::G4GAGTreeViewer
(G4VSceneHandler& sceneHandler, const G4String& name):
  G4VTreeViewer(sceneHandler, name) {
  // Make changes to view parameters for GAGTree...
  fVP.SetCulling(false);
  fDefaultVP.SetCulling(false);
}

G4GAGTreeViewer::~G4GAGTreeViewer() {}
