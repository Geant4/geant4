// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4RotationMatrix.hh,v 1.4 2000-11-24 09:50:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// ----------------------
// Class description:
//
// An "extended" rotation matrix class.
// The SetRotationMatrixByCol/Row() methods enables
// to define the rotation matrix by column/row vectors.
// The result matrix may be a matrix that does not
// represent a rotation transformation (!) as
// G3 "rotation" matrices can be a composition of
// rotation and reflection.  

// ----------------------

#ifndef G3TOG4ROTATION_HH
#define G3TOG4ROTATION_HH 1

#include "G4RotationMatrix.hh"
#include "globals.hh"

class G3toG4RotationMatrix : public G4RotationMatrix 
{

public:  // with description

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
