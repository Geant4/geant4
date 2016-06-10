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
// $Id: G3toG4RotationMatrix.cc 67982 2013-03-13 10:36:03Z gcosmo $

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
G3toG4RotationMatrix::SetRotationMatrixByCol(const G4ThreeVector& cl1,
                                             const G4ThreeVector& cl2,
                                             const G4ThreeVector& cl3)
{
  rxx = cl1.x();
  ryx = cl1.y();
  rzx = cl1.z();
  
  rxy = cl2.x();
  ryy = cl2.y();
  rzy = cl2.z();
  
  rxz = cl3.x();
  ryz = cl3.y();
  rzz = cl3.z();
  
}

void 
G3toG4RotationMatrix::SetRotationMatrixByRow(const G4ThreeVector& rw1,
                                             const G4ThreeVector& rw2,
                                             const G4ThreeVector& rw3)
{
  rxx = rw1.x();
  rxy = rw1.y();
  rxz = rw1.z();
  
  ryx = rw2.x();
  ryy = rw2.y();
  ryz = rw2.z();
  
  rzx = rw3.x();
  rzy = rw3.y();
  rzz = rw3.z();
  
}



