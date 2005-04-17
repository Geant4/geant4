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
// $Id: G4OpenGLFontBaseStore.cc,v 1.1 2005-04-17 16:08:43 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "G4OpenGLFontBaseStore.hh"

std::multimap<G4VViewer*,G4OpenGLFontBaseStore::FontInfo>
G4OpenGLFontBaseStore::fFontBaseMap;

void G4OpenGLFontBaseStore::AddFontBase
(G4VViewer* viewer,G4int fontBase,
 G4double size, const G4String& fontName) {
  fFontBaseMap.insert
    (std::make_pair(viewer, FontInfo(fontName, size, fontBase)));
}

G4int G4OpenGLFontBaseStore::GetFontBase(G4VViewer* viewer, G4double size) {
  G4int fontBase = -1;
  std::multimap<G4VViewer*,FontInfo>::const_iterator i, j;
  i = fFontBaseMap.find(viewer);
  if (i != fFontBaseMap.end()) {
    G4double sizeDiscrepancy = 9999.;
    while (i->first == viewer){
      G4double d = std::abs(size - i->second.fSize);
      if (d < sizeDiscrepancy) {
	sizeDiscrepancy = d;
	j = i;
      }
      ++i;
    }
    fontBase = j->second.fFontBase;
  }
  return fontBase;
}
