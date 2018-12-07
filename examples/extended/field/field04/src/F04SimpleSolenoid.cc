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
/// \file field/field04/src/F04SimpleSolenoid.cc
/// \brief Implementation of the F04SimpleSolenoid class
//

#include "globals.hh"

#include "G4GeometryManager.hh"

#include "F04GlobalField.hh"

#include "F04SimpleSolenoid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04SimpleSolenoid::F04SimpleSolenoid(G4double Bz, G4double fz,
                               G4LogicalVolume* lv,
                               G4ThreeVector c) : F04ElementField(c,lv)
{
  fBfield  = Bz;
  fFringeZ = fz;

  fFieldLength = 2.*((G4Tubs*)fVolume->GetSolid())->GetZHalfLength()+fFringeZ;
  fFieldRadius = ((G4Tubs*)fVolume->GetSolid())->GetOuterRadius();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04SimpleSolenoid::AddFieldValue(const G4double point[4],
                                            G4double field[6]) const
{
   G4ThreeVector global(point[0],point[1],point[2]);
   G4ThreeVector local;

   local = fGlobal2local.TransformPoint(global);

   if (IsOutside(local)) return;

   G4ThreeVector B(0.0,0.0,fBfield);

   B = fGlobal2local.Inverse().TransformAxis(B);

   field[0] += B[0];
   field[1] += B[1];
   field[2] += B[2];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool F04SimpleSolenoid::IsOutside(G4ThreeVector& local) const
{
//  EInside inside = tubs->Inside(local);
//  return (inside == kOutside);
  G4double r = std::sqrt(local.x()*local.x()+local.y()*local.y());
  return (r > fFieldRadius || std::fabs(local.z()) > fFieldLength/2.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool F04SimpleSolenoid::IsWithin(G4ThreeVector& local) const
{
//  EInside inside = tubs->Inside(local);
//  return (inside == kInside);
  G4double r = std::sqrt(local.x()*local.x()+local.y()*local.y());
  return (r < fFieldRadius && std::fabs(local.z()) < fFieldLength/2.0);
}
