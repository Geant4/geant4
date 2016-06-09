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
// $Id: G4OpenGLViewerDataStore.cc,v 1.1 2005/09/29 14:26:18 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4OpenGLViewerDataStore.hh"

std::map<G4VViewer*,G4bool> G4OpenGLViewerDataStore::fTransparencyEnabledMap;

void G4OpenGLViewerDataStore::SetTransparencyEnabled
(G4VViewer* viewer, G4bool transparency_enabled) {
  fTransparencyEnabledMap[viewer] = transparency_enabled;
}

G4bool G4OpenGLViewerDataStore::GetTransparencyEnabled(G4VViewer* viewer) {
  std::map<G4VViewer*,G4bool>::const_iterator i;
  i = fTransparencyEnabledMap.find(viewer);
  if (i != fTransparencyEnabledMap.end()) {
    return i->second;
  } else {
    return true;  // Default.
  }
}
