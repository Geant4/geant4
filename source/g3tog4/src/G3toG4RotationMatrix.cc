// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G3toG4RotationMatrix.cc,v 1.2 1999-12-05 17:50:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G3toG4RotationMatrix.hh"

G3toG4RotationMatrix::G3toG4RotationMatrix()
{
  rxx = 1;
  ryx = 0;
  rzx = 0;
  rxy = 0;
  ryy = 1;
  rzy = 0;
  rxz = 0;
  ryz = 0;
  rzz = 1;
}

G3toG4RotationMatrix::~G3toG4RotationMatrix()
{
    ;
}

void 
G3toG4RotationMatrix::SetRotationMatrixByCol(const G4ThreeVector& col1,
                                             const G4ThreeVector& col2,
                                             const G4ThreeVector& col3)
{
  rxx = col1.x();
  ryx = col1.y();
  rzx = col1.z();
  
  rxy = col2.x();
  ryy = col2.y();
  rzy = col2.z();
  
  rxz = col3.x();
  ryz = col3.y();
  rzz = col3.z();
  
}

void 
G3toG4RotationMatrix::SetRotationMatrixByRow(const G4ThreeVector& row1,
                                             const G4ThreeVector& row2,
                                             const G4ThreeVector& row3)
{
  rxx = row1.x();
  rxy = row1.y();
  rxz = row1.z();
  
  ryx = row2.x();
  ryy = row2.y();
  ryz = row2.z();
  
  rzx = row3.x();
  rzy = row3.y();
  rzz = row3.z();
  
}



