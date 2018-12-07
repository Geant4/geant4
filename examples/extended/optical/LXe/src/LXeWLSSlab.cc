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
/// \file optical/LXe/src/LXeWLSSlab.cc
/// \brief Implementation of the LXeWLSSlab class
//
//
#include "LXeWLSSlab.hh"
#include "LXeWLSFiber.hh"

#include "globals.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4SystemOfUnits.hh"

G4LogicalVolume* LXeWLSSlab::fScintSlab_log = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

LXeWLSSlab::LXeWLSSlab(G4RotationMatrix *pRot,
                       const G4ThreeVector &tlate,
                       G4LogicalVolume *pMotherLogical,
                       G4bool pMany,
                       G4int pCopyNo,
                       LXeDetectorConstruction* c)
  :G4PVPlacement(pRot,tlate,
                 new G4LogicalVolume(new G4Box("temp",1,1,1),
                                     G4Material::GetMaterial("Vacuum"),
                                     "temp",0,0,0),
                 "Slab",pMotherLogical,pMany,pCopyNo),fConstructor(c)
{
  CopyValues();
 
  G4double slab_x = fScint_x/2.;
  G4double slab_y = fScint_y/2.;
 
  G4Box* ScintSlab_box = new G4Box("Slab",slab_x,slab_y,fSlab_z);
 
  fScintSlab_log
      = new G4LogicalVolume(ScintSlab_box,
                            G4Material::GetMaterial("Polystyrene"),
                            "Slab",0,0,0);
 
  G4double spacing = 2*slab_y/fNfibers;
 
  G4RotationMatrix* rm = new G4RotationMatrix();
  rm->rotateY(90*deg);
 
  //Place fibers
  for(G4int i=0;i<fNfibers;i++){
     G4double Y=-(spacing)*(fNfibers-1)*0.5 + i*spacing;
     new LXeWLSFiber(rm,G4ThreeVector(0.,Y,0.),fScintSlab_log,false,0,
                     fConstructor);
  }
 
  SetLogicalVolume(fScintSlab_log);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void LXeWLSSlab::CopyValues(){
 
  fScint_x=fConstructor->GetScintX();
  fScint_y=fConstructor->GetScintY();
  fScint_z=fConstructor->GetScintZ();
  fNfibers=fConstructor->GetNFibers();
  fSlab_z=fConstructor->GetSlabZ();
}
