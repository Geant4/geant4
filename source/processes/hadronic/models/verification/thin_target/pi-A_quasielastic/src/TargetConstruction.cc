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
// $Id: TargetConstruction.cc,v 1.1 2003-07-31 01:22:03 dwright Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
#include "TargetConstruction.hh"
#include "TargetMessenger.hh"

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"


TargetConstruction::TargetConstruction()
 :WorldMaterial(0), solidWorld(0), logicWorld(0), physWorld(0),
  TargetMaterial(0), solidTarg(0), logicTarg(0), physTarg(0)
{
  // Default target parameters
  TargetThickness = 1.674*mm;
  TargetRadius = 50.0*mm;
  NumRequestedEvents = 1;

  // create commands for interactive definition of the calorimeter  
  tgtMessenger = new TargetMessenger(this);
}


TargetConstruction::~TargetConstruction()
{ 
  delete tgtMessenger;
}


G4VPhysicalVolume* TargetConstruction::Construct()
{
  DefineMaterials();
  return ConstructTarget();
}


void TargetConstruction::DefineMaterials()
{ 
  new G4Material("Deuterium", 1.,  2.01*g/mole, 0.162*g/cm3);
  new G4Material("Lithium6",  3.,  6.00*g/mole, 0.457*g/cm3);
  new G4Material("Beryllium", 4.,  9.01*g/mole, 1.85*g/cm3);
  new G4Material("Carbon",    6., 12.01*g/mole, 2.26*g/cm3);
  new G4Material("Aluminum", 13., 26.98*g/mole, 2.70*g/cm3);
  new G4Material("Calcium",  20., 40.08*g/mole, 1.55*g/cm3);
  new G4Material("Tin",      50.,118.69*g/mole, 8.82*g/cm3);
  new G4Material("Copper",   29., 63.55*g/mole, 8.96*g/cm3);
  new G4Material("Zirconium",40., 91.22*g/mole, 6.51*g/cm3);
  new G4Material("Tungsten", 74.,183.85*g/mole,19.30*g/cm3);
  new G4Material("Platinum", 78.,195.08*g/mole,21.45*g/cm3);
  new G4Material("Lead",     82.,207.19*g/mole,11.35*g/cm3);
  G4Material* U = new G4Material("Uranium", 92.,238.03*g/mole,18.95*g/cm3);

  G4Material* Vacuum = new G4Material("Galactic", 1., 1.01*g/mole,
           universe_mean_density,kStateGas,2.73*kelvin,3.e-18*pascal);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //default materials of the target
  TargetMaterial = U;
  WorldMaterial  = Vacuum;
}

  
G4VPhysicalVolume* TargetConstruction::ConstructTarget()
{
  //     
  // World
  //
  WorldLength = TargetThickness + 1.0*cm;
  G4double WorldRadius = TargetRadius + 0.5*cm;
  G4cout << " TargetCon: WorldLength = " << WorldLength/mm << G4endl;
  
  solidWorld = new G4Tubs("World", 0., WorldRadius, WorldLength/2, 
                                                            0., 360.0*deg); 
  logicWorld = new G4LogicalVolume(solidWorld, WorldMaterial, "World");
  physWorld = new G4PVPlacement(0, 0, "World", logicWorld, 0, false, 0);
  
  //                               
  // Target
  //  

  solidTarg = new G4Tubs("Target", 0., TargetRadius, TargetThickness/2., 
                                                            0., 360.0*deg);
  logicTarg = new G4LogicalVolume(solidTarg, TargetMaterial,"Target");
  physTarg = new G4PVPlacement(0,0,"Target",logicTarg,physWorld,false,0);
  
  //                                        
  // Visualization attributes
  //
  logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleVisAtt->SetVisibility(true);
  logicTarg->SetVisAttributes(simpleVisAtt);
  
  //
  // Always return the physical World
  //
  PrintTargParameters();
  return physWorld;
}


void TargetConstruction::PrintTargParameters()
{
  G4cout << "\n-----------------------------------------------------------"
         << "\n---> Target is " << TargetMaterial->GetName() 
         << " " << TargetThickness/mm << "mm thick and radius "  
         << TargetRadius/mm << " mm "
         << "\n-----------------------------------------------------------\n";
}


void TargetConstruction::SetTargetMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialChoice);     
  if(pttoMaterial) {
    TargetMaterial = pttoMaterial;
    logicTarg->SetMaterial(pttoMaterial); 
    PrintTargParameters();
  } else {
    G4cout << " mat ptr 0 " << G4endl;
  }
}


void TargetConstruction::SetTargetThickness(G4double val)
{
  TargetThickness = val;
  PrintTargParameters();
}  


void TargetConstruction::SetTargetRadius(G4double val)
{
  TargetRadius = val;
  PrintTargParameters();
}  


void TargetConstruction::SetNumRequestedEvents(G4int val)
{
  NumRequestedEvents = val;
  PrintTargParameters();
}  


void TargetConstruction::SetMagField(G4double fieldValue)
{
}

  
void TargetConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructTarget());
}


