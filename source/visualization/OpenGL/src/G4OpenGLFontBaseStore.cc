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
// $Id: G4OpenGLFontBaseStore.cc,v 1.2 2005/05/27 10:56:06 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//

#include "G4OpenGLFontBaseStore.hh"

std::map<G4VViewer*,std::vector<G4OpenGLFontBaseStore::FontInfo> >
G4OpenGLFontBaseStore::fFontBaseMap;

void G4OpenGLFontBaseStore::AddFontBase
(G4VViewer* viewer,G4int fontBase,
 G4double size, const G4String& fontName) {
  fFontBaseMap[viewer].push_back(FontInfo(fontName, size, fontBase));
}

G4int G4OpenGLFontBaseStore::GetFontBase(G4VViewer* viewer, G4double size) {
  G4int fontBase = -1;
  std::map<G4VViewer*,std::vector<FontInfo> >::const_iterator i;
  i = fFontBaseMap.find(viewer);
  if (i != fFontBaseMap.end()) {
    G4double sizeDiscrepancy = 9999.;
    std::vector<FontInfo>::const_iterator j, k;
    for (j = i->second.begin(); j != i->second.end(); ++j) {
      G4double d = std::abs(size - j->fSize);
      if (d < sizeDiscrepancy) {
	sizeDiscrepancy = d;
	k = j;
      }
    }
    fontBase = k->fFontBase;
  }
  return fontBase;
}
