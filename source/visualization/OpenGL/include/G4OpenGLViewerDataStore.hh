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
// $Id: G4OpenGLViewerDataStore.hh,v 1.1 2005/09/29 14:26:17 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// J.Allison  Apr 2005.

// Class Description:
// Facilitates sharing of fontdat between scene handler and viewers.

#ifndef G4OPENGLVIEWERDATASTORE_HH
#define G4OPENGLVIEWERDATASTORE_HH

#include "globals.hh"
#include <map>
#include <vector>

class G4VViewer;

class G4OpenGLViewerDataStore {
public:
  static void SetTransparencyEnabled(G4VViewer*, G4bool transparency_enabled);
  static G4bool GetTransparencyEnabled(G4VViewer*);
private:
  static std::map<G4VViewer*,G4bool> fTransparencyEnabledMap;
};

#endif
