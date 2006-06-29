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
// $Id: G3toG4RotationMatrix.cc,v 1.4 2006-06-29 18:13:31 gunter Exp $
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



