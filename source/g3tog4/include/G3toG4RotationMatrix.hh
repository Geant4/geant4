//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G3toG4RotationMatrix.hh,v 1.5 2001-07-11 09:58:58 gunter Exp $
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
