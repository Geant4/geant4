// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3CalcParams.hh,v 1.2 1999-12-05 17:50:00 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G3CalcParams
//
// Class inheriting from G4VPVParameterisation to provide GSDVN
// functionality
//
// History:
// 03.09.96 T. Wenaus  Initial version

#ifndef G3CALCPARAMS_HH
#define G3CALCPARAMS_HH

#include <math.h>
#include "globals.hh"
#include "G4VPVParameterisation.hh"

class G3CalcParams: public G4VPVParameterisation
{
public:
  G3CalcParams(G4int ndiv, G4double width, G4double offset)
    { fNdiv = ndiv; fWidth = width; fOffset = offset;}

  void ComputeTransformation(const G4int n,
                             G4VPhysicalVolume *pRep) const;
  
  void ComputeDimensions(G4Box &pBox,
                         const G4int,
                         const G4VPhysicalVolume *pRep) const;
  
  void ComputeDimensions(G4Tubs &,
                         const G4int,
                         const G4VPhysicalVolume *) const;
                         
  void ComputeDimensions(G4Trd &,
                         const G4int,
                         const G4VPhysicalVolume *) const;
 
  void ComputeDimensions(G4Trap &,
                         const G4int,
                         const G4VPhysicalVolume *) const; 
 
  void ComputeDimensions(G4Cons &,
                         const G4int,
                         const G4VPhysicalVolume *) const;
 
  void ComputeDimensions(G4Sphere &,
                         const G4int,
                         const G4VPhysicalVolume *) const;

  void ComputeDimensions(G4Torus &,
                         const G4int,
                         const G4VPhysicalVolume *) const;

  void ComputeDimensions(G4Para &,
                         const G4int,
                         const G4VPhysicalVolume *) const;
  
  void ComputeDimensions(G4Hype &,
                         const G4int,
                         const G4VPhysicalVolume *) const;

private:
  G4int fNdiv;
  G4double fWidth;
  G4double fOffset;
  
};
#endif
