// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParametrizedBox.hh,v 1.2 1999-12-15 14:54:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Original class written by Hans-Peter Wellisch.

#include "G4VPVParameterisation.hh"
 class G4ParametrizedBox: public G4VPVParameterisation
 {
  virtual void ComputeTransformation(const G4int n,
                                     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0,(n-1)*15*m,0));
  }

  virtual void ComputeDimensions(G4Box &pBox,
                                 const G4int n,
                                 const G4VPhysicalVolume* pRep) const
  {
    pBox.SetXHalfLength(10*m*n);
    pBox.SetYHalfLength(5*m*n);
    pBox.SetZHalfLength(5*m*n);
  }
 };
