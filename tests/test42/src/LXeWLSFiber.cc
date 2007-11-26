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
#include "LXeWLSFiber.hh"
#include "globals.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

G4LogicalVolume* LXeWLSFiber::clad2_log=NULL;

LXeWLSFiber::LXeWLSFiber(G4RotationMatrix *pRot,
			     const G4ThreeVector &tlate,
			     G4LogicalVolume *pMotherLogical,
			     G4bool pMany,
			     G4int pCopyNo,
			     LXeDetectorConstruction* c)
  :G4PVPlacement(pRot,tlate,
		 new G4LogicalVolume(new G4Box("temp",1,1,1),
				     G4Material::GetMaterial("Vacuum"),
				     "temp",0,0,0),
		 "Cladding2",pMotherLogical,pMany,pCopyNo),constructor(c)
{
  CopyValues();
  
  if(!clad2_log || updated){
    // The Fiber
    //
    G4Tubs* Fiber_tube = 
      new G4Tubs("Fiber",fiber_rmin,fiber_rmax,fiber_z,fiber_sphi,fiber_ephi);
    
    G4LogicalVolume* Fiber_log = 
      new G4LogicalVolume(Fiber_tube,G4Material::GetMaterial("PMMA"),
			  "Fiber",0,0,0);
    
    // Cladding (first layer)
    //
    G4Tubs* clad1_tube = 
      new G4Tubs("Cladding1",clad1_rmin,clad1_rmax,clad1_z,clad1_sphi,
		 clad1_ephi);
    
    G4LogicalVolume* clad1_log = 
      new G4LogicalVolume(clad1_tube,G4Material::GetMaterial("Pethylene"),
			  "Cladding1",0,0,0);
     
    // Cladding (second layer)
    // 
    G4Tubs* clad2_tube = 
      new G4Tubs("Cladding2",clad2_rmin,clad2_rmax,clad2_z,clad2_sphi,
		 clad2_ephi);
    
    clad2_log = 
      new G4LogicalVolume(clad2_tube,G4Material::GetMaterial("fPethylene"),
			  "Cladding2",0,0,0);
    
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Fiber_log,
		      "Fiber", clad1_log,false,0);
    new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),clad1_log,
		      "Cladding1",clad2_log,false,0);
  }
  
  SetLogicalVolume(clad2_log);
}

void LXeWLSFiber::CopyValues(){
  updated=constructor->GetUpdated();

  fiber_rmin = 0.00*cm;    
  fiber_rmax = 0.10*cm;    
  fiber_z    = constructor->GetScintX()/2;
  fiber_sphi = 0.00*deg;
  fiber_ephi = 360.*deg;
  
  clad1_rmin = 0.;// fiber_rmax;
    clad1_rmax = fiber_rmax + 0.015*fiber_rmax;    
  
  clad1_z    = fiber_z;
  clad1_sphi = fiber_sphi;
  clad1_ephi = fiber_ephi; 
  
  clad2_rmin = 0.;//clad1_rmax;
  clad2_rmax = clad1_rmax + 0.015*fiber_rmax;    

  clad2_z    = fiber_z;
  clad2_sphi = fiber_sphi;
  clad2_ephi = fiber_ephi;

}
