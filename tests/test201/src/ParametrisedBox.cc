// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ParametrisedBox.cc,v 1.1 1999-01-08 16:35:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Test of parametrisation.  Based on Hans-Peter Wellisch's test.

#include "ParametrisedBox.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"

void ParametrisedBox::ComputeTransformation
(const G4int n,
 G4VPhysicalVolume* pRep) const
{
  pRep->SetTranslation (G4ThreeVector (0, n * 2 * m, 0));
}

void ParametrisedBox::ComputeDimensions
(G4Box& box,
 const G4int n,
 const G4VPhysicalVolume* pRep) const
{
  box.SetXHalfLength ( 0.5 * m * (n + 1));
  box.SetYHalfLength ( 1 * m * (n + 1));
  box.SetZHalfLength ( 1.5 * m * (n + 1));
}
