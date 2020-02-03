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
/// \file runAndEvent/RE04/src/RE04ParallelWorldParam.cc
/// \brief Implementation of the RE04ParallelWorldParam class
//
//
#include "RE04ParallelWorldParam.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04ParallelWorldParam::RE04ParallelWorldParam()
 : G4VPVParameterisation(), 
   fWater(0), fPb(0)
{
  fWater = G4Material::GetMaterial("G4_WATER");
  fPb = G4Material::GetMaterial("G4_Pb");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04ParallelWorldParam::~RE04ParallelWorldParam()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04ParallelWorldParam::ComputeTransformation(const G4int copyNo,
                                     G4VPhysicalVolume* physVol) const
{
  G4double x = copyNo ? -10.0*cm : 10.0*cm;
  physVol->SetTranslation(G4ThreeVector(x,x,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4Material* RE04ParallelWorldParam::ComputeMaterial(const G4int copyNo,
      G4VPhysicalVolume* /*currentVol*/,const G4VTouchable* /*parentTouch*/)
{ return copyNo ? fWater: fPb; }

