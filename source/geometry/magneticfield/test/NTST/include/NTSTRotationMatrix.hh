// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTRotationMatrix.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
//
#ifndef _NTSTROTATION_
#define _NTSTROTATION_

#include "G4RotationMatrix.hh"
#include "globals.hh"

class NTSTRotationMatrix : public G4RotationMatrix 
{
public:
  NTSTRotationMatrix();

  void SetRotationMatrixByCol(const G4ThreeVector& Col1,
                              const G4ThreeVector& Col2,
                              const G4ThreeVector& Col3);
    
  void SetRotationMatrixByRow(const G4ThreeVector& Row1,
                              const G4ThreeVector& Row2,
                              const G4ThreeVector& Row3);
    
    ~NTSTRotationMatrix();
};

#endif    
