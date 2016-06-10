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
// $Id: G4OpenGLFontBaseStore.cc 66373 2012-12-18 09:41:34Z gcosmo $
//

#include "G4OpenGLFontBaseStore.hh"

std::map<G4VViewer*,std::vector<G4OpenGLFontBaseStore::FontInfo> >
G4OpenGLFontBaseStore::fFontBaseMap;

void G4OpenGLFontBaseStore::AddFontBase
(G4VViewer* viewer, G4int fontBase,
 G4double size, const G4String& fontName,
 G4int width) {
  fFontBaseMap[viewer].push_back(FontInfo(fontName, size, fontBase, width));
}

const G4OpenGLFontBaseStore::FontInfo&
G4OpenGLFontBaseStore::GetFontInfo
(G4VViewer* viewer, G4double size) {
  std::map<G4VViewer*,std::vector<FontInfo> >::const_iterator i =
    fFontBaseMap.find(viewer);
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
    return *k;
  } else {
    // No font found for requested viewer.
    static const FontInfo nullFontInfo;  // Default struct, fontBase = -1;
    return nullFontInfo;
  }
}
