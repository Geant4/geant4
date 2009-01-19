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
// $Id: G4OpenGLFontBaseStore.cc,v 1.4 2009-01-19 16:53:42 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
