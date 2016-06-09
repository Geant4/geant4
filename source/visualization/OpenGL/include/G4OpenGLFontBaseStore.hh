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
// $Id: G4OpenGLFontBaseStore.hh,v 1.2 2005/05/27 10:56:06 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// J.Allison  Apr 2005.

// Class Description:
// Facilitates sharing of font bases between scene handler and viewers.

#ifndef G4OPENGLFONTBASESTORE_HH
#define G4OPENGLFONTBASESTORE_HH

#include "globals.hh"
#include <map>
#include <vector>

class G4VViewer;

class G4OpenGLFontBaseStore {
public:
  static void AddFontBase(G4VViewer*, G4int fontBase,
			  G4double size, const G4String& fontName);
  static G4int GetFontBase(G4VViewer*, G4double size);
private:
  struct FontInfo {
    FontInfo():
      fFontName(""), fSize(0), fFontBase(-1) {}
    FontInfo(const G4String& fontName, G4double size, G4int fontBase):
      fFontName(fontName), fSize(size), fFontBase(fontBase) {}
    G4String fFontName;
    G4double fSize;  // In terms of G4VMarker Screen Size.
    G4int fFontBase;
  };
  static std::map<G4VViewer*,std::vector<FontInfo> > fFontBaseMap;
};

#endif
