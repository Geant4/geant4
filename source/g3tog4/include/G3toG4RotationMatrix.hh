// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4RotationMatrix.hh,v 1.2 1999-12-05 17:50:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef _G3toG4ROTATION_
#define _G3toG4ROTATION_

#include "G4RotationMatrix.hh"
#include "globals.hh"

class G3toG4RotationMatrix : public G4RotationMatrix 
{
public:
  G3toG4RotationMatrix();

  void SetRotationMatrixByCol(const G4ThreeVector& Col1,
                              const G4ThreeVector& Col2,
                              const G4ThreeVector& Col3);
    
  void SetRotationMatrixByRow(const G4ThreeVector& Row1,
                              const G4ThreeVector& Row2,
                              const G4ThreeVector& Row3);
    
    ~G3toG4RotationMatrix();
};

#endif    
