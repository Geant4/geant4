// $Id: A01CellParameterisation.cc,v 1.1 2002-11-13 07:22:43 duns Exp $
// --------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "A01CellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

A01CellParameterisation::A01CellParameterisation()
{
  for(int copyNo=0;copyNo<80;copyNo++)
  {
    G4int column = copyNo / 4;
    G4int row = copyNo % 4;
    xCell[copyNo] = (column-9)*15.*cm - 7.5*cm;
    yCell[copyNo] = (row-1)*15*cm - 7.5*cm;
  }
}

A01CellParameterisation::~A01CellParameterisation()
{;}

void A01CellParameterisation::ComputeTransformation
(const G4int copyNo,G4VPhysicalVolume *physVol) const
{
  physVol->SetTranslation(G4ThreeVector(xCell[copyNo],yCell[copyNo],0.));
}

