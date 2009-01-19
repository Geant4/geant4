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
// $Id: G4OpenGLBitMapStore.cc,v 1.5 2009-01-19 16:53:42 lgarnier Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  6th January 2007

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLBitMapStore.hh"

namespace G4OpenGLBitMapStore {

std::map<Key, GLubyte*> fStore;

void SetBit(G4int i, G4int j, GLubyte* bitmap, G4int byteColumns)
{
  G4int k = i / 8 + byteColumns * j;
  bitmap[k] |= 0x80 >> i % 8;
}

const GLubyte* GetBitMap(Shape shape, G4double& size, G4bool filled)
{
  // Rationalise size
  G4int bitSize = G4int(std::abs(size) + 0.49);
  if (bitSize < 1) bitSize = 1;

  // Modify size in calling function (passed by reference).
  size = bitSize;

  Key key(shape, bitSize, filled);

  // If previously encountered, return.
  if (fStore.find(key) != fStore.end()) return fStore[key];

  // Create and store new bitmap.
  G4int byteColumns = (G4int(bitSize - 1) / 8) + 1;
  G4int byteRows = bitSize;
  G4int nBytes = byteColumns*byteRows;
  GLubyte* bitmap = new GLubyte[nBytes];
  fStore[key] = bitmap;

  if (shape == square && filled) {  // Set all bits and return.
    for (G4int i = 0; i < nBytes; ++i) bitmap[i] = 0xff;
    return bitmap;
  } else {  // Clear ready for setting bits.
    for (G4int i = 0; i < nBytes; ++i) bitmap[i] = 0x00;
  }

  // Fill bitmap.
  GLint linewidth;
  glGetIntegerv(GL_LINE_WIDTH, &linewidth);
  if (linewidth > bitSize) linewidth = bitSize;
  G4double outer = bitSize / 2.;
  G4double outer2 = outer * outer;
  G4double inner = bitSize / 2. - linewidth;
  G4double inner2 = inner * inner;
  for (G4int i = 0; i < bitSize; ++i) {
    for (G4int j = 0; j < bitSize; ++j) {
      G4double x = i - (bitSize - 1) * 0.5;
      G4double y = j - (bitSize - 1) * 0.5;
      G4double r2 = x * x + y * y;
      if (shape == circle) {
	if (filled) {
	  if (r2 < outer2) SetBit(i, j, bitmap, byteColumns);
	} else {
	  if (r2 < outer2 && r2 > inner2) SetBit(i, j, bitmap, byteColumns);
	}
      } else {  // Unfilled square.
	if (x < -inner || x > inner || y < -inner || y > inner)
	  SetBit(i, j, bitmap, byteColumns);
      }
    }
  }

  return bitmap;
}

} // End namespace G4OpenGLBitMapStore.

#endif
