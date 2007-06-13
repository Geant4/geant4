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
// $Id: B02PSScoringDetectorConstruction.cc,v 1.1 2007-06-13 13:31:41 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include <sstream>

#include "B02PSScoringDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"



B02PSScoringDetectorConstruction::B02PSScoringDetectorConstruction(G4String worldName)
  :G4VUserParallelWorld(worldName),fLogicalVolumeVector() //ASO
{
  //  Construct();
}

B02PSScoringDetectorConstruction::~B02PSScoringDetectorConstruction()
{
  fLogicalVolumeVector.clear(); //ASO
}

void B02PSScoringDetectorConstruction::Construct()
{  

  G4cout << " constructing parallel world " << G4endl;

  //GetWorld methods create a clone of the mass world to the parallel world (!)
  // via the transportation manager
  ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  fLogicalVolumeVector.push_back(worldLogical); //ASO

  G4String name("none");
  G4double density(universe_mean_density), temperature(0), pressure(0);

  name = "Galactic";
  density     = universe_mean_density;            //from PhysicalConstants.h
  pressure    = 3.e-18*pascal;
  temperature = 2.73*kelvin;
  G4cout << density << " " << kStateGas << G4endl;
  G4Material *Galactic = 
    new G4Material(name, 1., 1.01*g/mole, density,
                   kStateGas,temperature,pressure);


  //////////////////////////////////
  // parallel world cylinder volume
  //////////////////////////////////

  // parallel world solid larger than in the mass geometry

  //G4double innerRadiusCylinder = 0*cm;
  //G4double outerRadiusCylinder = 110*cm;
  //G4double hightCylinder       = 110*cm;
  //G4double startAngleCylinder  = 0*deg;
  //G4double spanningAngleCylinder    = 360*deg;

  /*
  G4Tubs *worldCylinder = new G4Tubs("worldCylinder",
                                     innerRadiusCylinder,
                                     outerRadiusCylinder,
                                     hightCylinder,
                                     startAngleCylinder,
                                     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, Galactic, "worldCylinder_log");

  name = "parallelWorld";
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		  name, 0, false, 0);

  */
  //fPVolumeStore.AddPVolume(G4GeometryCell(*ghostWorld, 0));




  // creating 18 slobs of 10 cm thicknes

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 100*cm;
  G4double hightShield       = 5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               hightShield,
                               startAngleShield,
                               spanningAngleShield);
  
  // logical parallel cells

  G4LogicalVolume *aShield_log = 
    new G4LogicalVolume(aShield, Galactic, "aShield_log");
  fLogicalVolumeVector.push_back(aShield_log); //ASO

  // physical parallel cells

  G4int i = 1;
  G4double startz = -85*cm; 
  for (i=1; i<=18; ++i) {
   
    name = GetCellName(i);
    
    G4double pos_x = 0*cm;
    G4double pos_y = 0*cm;
    G4double pos_z = startz + (i-1) * (2*hightShield);
    G4VPhysicalVolume *pvol = 
      new G4PVPlacement(0, 
			G4ThreeVector(pos_x, pos_y, pos_z),
			aShield_log, 
			name, 
			worldLogical, 
			false, 
			i); //ASO
    //G4GeometryCell cell(*pvol, i);  //ASO
    //fPVolumeStore.AddPVolume(cell);
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*cm;
  outerRadiusShield = 100*cm;  //ASO
  hightShield       = 5*cm;    // ASO 
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;

  G4Tubs *aRest = new G4Tubs("Rest",
			     innerRadiusShield,
			     outerRadiusShield,
			     hightShield,
			     startAngleShield,
			     spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, Galactic, "aRest_log");
  name = GetCellName(19);
  fLogicalVolumeVector.push_back(aRest_log); //ASO
    
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  G4double pos_z = 95*cm;
  G4VPhysicalVolume *pvol = 
    new G4PVPlacement(0, 
		      G4ThreeVector(pos_x, pos_y, pos_z),
		      aRest_log, 
		      name, 
		      worldLogical, 
		      false, 
		      19); // i = 19  ASO
  //G4GeometryCell cell(*pvol, i); // ASO
  //fPVolumeStore.AddPVolume(cell);

  SetSensitive(); //ASO

}

//const G4VPhysicalVolume &B02PSScoringDetectorConstruction::
//GetPhysicalVolumeByName(const G4String& name) const {
//  return *fPVolumeStore.GetPVolume(name);
//}


//G4String B02PSScoringDetectorConstruction::ListPhysNamesAsG4String(){
//  G4String names(fPVolumeStore.GetPNames());
//  return names;
//}


G4String B02PSScoringDetectorConstruction::GetCellName(G4int i) {
  std::ostringstream os;
  os << "cell_";
  if (i<10) {
    os << "0";
  }
  os << i;
  G4String name = os.str();
  return name;
}

//G4GeometryCell B02PSScoringDetectorConstruction::GetGeometryCell(G4int i){
//  G4String name(GetCellName(i));
//  const G4VPhysicalVolume *p=0;
//  p = fPVolumeStore.GetPVolume(name);
//  if (p) {
//    G4cout << " returning GetGeometryCell " << G4endl;
//    return G4GeometryCell(*p,i); //ASO
//  }
//  else {
//    G4cout << "B02PSScoringDetectorConstruction::GetGeometryCell: couldn't get G4GeometryCell" << G4endl;
//    return G4GeometryCell(*ghostWorld,-2);
//  }
//}


// G4VPhysicalVolume *B02PSScoringDetectorConstruction::GetWorldVolume() const{
//   return *fWorldVolume;
// }


G4VPhysicalVolume *B02PSScoringDetectorConstruction::GetWorldVolume() {
   return ghostWorld;
}


G4VPhysicalVolume &B02PSScoringDetectorConstruction::GetWorldVolumeAddress() const{
  return *ghostWorld;
}

void B02PSScoringDetectorConstruction::SetSensitive(){

  //  -------------------------------------------------
  //   The collection names of defined Primitives are
  //   0       PhantomSD/Collisions
  //   1       PhantomSD/CollWeight
  //   2       PhantomSD/Population
  //   3       PhantomSD/TrackEnter
  //   4       PhantomSD/SL
  //   5       PhantomSD/SLW
  //   6       PhantomSD/SLWE
  //   7       PhantomSD/SLW_V
  //   8       PhantomSD/SLWE_V
  //  -------------------------------------------------


  //================================================
  // Sensitive detectors : MultiFunctionalDetector
  //================================================
  //
  //  Sensitive Detector Manager.
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //
  // Sensitive Detector Name
  G4String phantomSDname = "PhantomSD";

  //------------------------
  // MultiFunctionalDetector
  //------------------------
  //
  // Define MultiFunctionalDetector with name.
  G4MultiFunctionalDetector* MFDet = new G4MultiFunctionalDetector(phantomSDname);
  SDman->AddNewDetector( MFDet );                 // Register SD to SDManager


  G4String fltName,particleName;
  G4SDParticleFilter* neutronFilter = 
      new G4SDParticleFilter(fltName="neutronFilter", particleName="neutron");

  MFDet->SetFilter(neutronFilter);


  for (std::vector<G4LogicalVolume *>::iterator it =  fLogicalVolumeVector.begin();
       it != fLogicalVolumeVector.end(); it++){
      (*it)->SetSensitiveDetector(MFDet);
  }

  G4String psName;
  G4PSNofCollision*   scorer0 = new G4PSNofCollision(psName="Collisions");  
  MFDet->RegisterPrimitive(scorer0);


  G4PSNofCollision*   scorer1 = new G4PSNofCollision(psName="CollWeight");  
  scorer1->Weighted(true);
  MFDet->RegisterPrimitive(scorer1);


  G4PSPopulation*   scorer2 = new G4PSPopulation(psName="Population");  
  MFDet->RegisterPrimitive(scorer2);

  G4PSTrackCounter* scorer3 = new G4PSTrackCounter(psName="TrackEnter",fCurrent_In);  
  MFDet->RegisterPrimitive(scorer3);

  G4PSTrackLength* scorer4 = new G4PSTrackLength(psName="SL");  
  MFDet->RegisterPrimitive(scorer4);

  G4PSTrackLength* scorer5 = new G4PSTrackLength(psName="SLW");  
  scorer5->Weighted(true);
  MFDet->RegisterPrimitive(scorer5);

  G4PSTrackLength* scorer6 = new G4PSTrackLength(psName="SLWE");  
  scorer6->Weighted(true);
  scorer6->MultiplyKineticEnergy(true);
  MFDet->RegisterPrimitive(scorer6);

  G4PSTrackLength* scorer7 = new G4PSTrackLength(psName="SLW_V");  
  scorer7->Weighted(true);
  scorer7->DivideByVelocity(true);
  MFDet->RegisterPrimitive(scorer7);

  G4PSTrackLength* scorer8 = new G4PSTrackLength(psName="SLWE_V");  
  scorer8->Weighted(true);
  scorer8->MultiplyKineticEnergy(true);
  scorer8->DivideByVelocity(true);
  MFDet->RegisterPrimitive(scorer8);

}
