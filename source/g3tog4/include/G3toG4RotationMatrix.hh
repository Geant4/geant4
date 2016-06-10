//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G3toG4RotationMatrix.hh 67982 2013-03-13 10:36:03Z gcosmo $
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
