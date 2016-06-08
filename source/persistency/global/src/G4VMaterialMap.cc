// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VMaterialMap.cc,v 1.1 2000/06/09 12:53:51 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

#include "G4VMaterialMap.hh"

G4VMaterialMap* G4VMaterialMap::fMaterialMap = 0;

G4VMaterialMap* G4VMaterialMap::GetMaterialMap()
{
  return fMaterialMap;
}

G4VMaterialMap::G4VMaterialMap()
{
  fMaterialMap = this;
}

G4VMaterialMap::~G4VMaterialMap()
{
  fMaterialMap = NULL;
}


