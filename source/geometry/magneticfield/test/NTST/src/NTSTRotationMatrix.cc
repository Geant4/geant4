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
#include "NTSTRotationMatrix.hh"

NTSTRotationMatrix::NTSTRotationMatrix()
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

NTSTRotationMatrix::~NTSTRotationMatrix()
{
    ;
}

void 
NTSTRotationMatrix::SetRotationMatrixByCol(const G4ThreeVector& c1,
					   const G4ThreeVector& c2,
					   const G4ThreeVector& c3)
{
  rxx = c1.x();
  ryx = c1.y();
  rzx = c1.z();
  
  rxy = c2.x();
  ryy = c2.y();
  rzy = c2.z();
  
  rxz = c3.x();
  ryz = c3.y();
  rzz = c3.z();
}

void 
NTSTRotationMatrix::SetRotationMatrixByRow(const G4ThreeVector& r1,
					   const G4ThreeVector& r2,
					   const G4ThreeVector& r3)
{
  rxx = r1.x();
  rxy = r1.y();
  rxz = r1.z();
  
  ryx = r2.x();
  ryy = r2.y();
  ryz = r2.z();
  
  rzx = r3.x();
  rzy = r3.y();
  rzz = r3.z();
}



