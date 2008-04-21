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
// $Id: Tst33ParallelGeometry.cc,v 1.15 2008-04-21 09:00:03 ahoward Exp $
// GEANT4 tag 
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// Tst33ParallelGeometry.cc
//
// ----------------------------------------------------------------------

#include "Tst33ParallelGeometry.hh"

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GeometryCell.hh"

//ASO
// For Primitive Scorers
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDParticleFilter.hh"
#include "G4PSNofCollision.hh"
//#include "G4PSPopulation.hh"
#include "G4PSTrackCounter.hh"
#include "G4PSTrackLength.hh"


Tst33ParallelGeometry::Tst33ParallelGeometry(G4String worldName, G4VPhysicalVolume* ghostworld)
  :
  G4VUserParallelWorld(worldName),
  fLogicalVolumeVector(),
  fGalactic(0)
{
  worldVolumeName = worldName;
  ghostWorld = ghostworld;
  //  Construct();
}
Tst33ParallelGeometry::~Tst33ParallelGeometry(){}

void Tst33ParallelGeometry::Construct(){

  G4cout << " constructing parallel world " << G4endl;

  //GetWorld methods create a clone of the mass world to the parallel world (!)
  // via the transportation manager
  //x  
  ghostWorld = GetWorld();
  //  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();
  //  fLogicalVolumeVector.push_back(worldLogical); //ASO
  G4LogicalVolume* worldCylinder_log = ghostWorld->GetLogicalVolume();
  fLogicalVolumeVector.push_back(worldCylinder_log); //ASO


  fGalactic = fMaterialFactory.CreateGalactic();

  //////////////////////////////////
  // parallel world cylinder volume
  //////////////////////////////////

  // parallel world solid larger than in the mass geometry

  // don't create parallel world as it's now synonymous with physical world
  /*
  G4double innerRadiusCylinder = 0*cm;
  G4double outerRadiusCylinder = scaling*110*cm;
  G4double halfheightCylinder       = scaling*110*cm;
  G4double startAngleCylinder  = 0*deg;
  G4double spanningAngleCylinder    = 360*deg;

  G4Tubs *worldCylinder = new G4Tubs(worldVolumeName,
                                     innerRadiusCylinder,
                                     outerRadiusCylinder,
                                     halfheightCylinder,
                                     startAngleCylinder,
                                     spanningAngleCylinder);

  // logical world

  G4LogicalVolume *worldCylinder_log = 
    new G4LogicalVolume(worldCylinder, fGalactic, "worldCylinder_log");

  G4String name("parallelWorld");
  fWorldVolume = new 
    G4PVPlacement(0, G4ThreeVector(0,0,0), worldCylinder_log,
		  name, 0, false, 0);
  if (!fWorldVolume) {
    G4Exception("Tst33ParallelGeometry::Construct(): new failed to create G4PVPlacement!");
  }
  */

  //  fPVolumeStore.AddPVolume(G4GeometryCell(*fWorldVolume, -1));
  //xxx  fPVolumeStore.AddPVolume(G4GeometryCell(*ghostWorld, -1));




  // creating 18 slobs of 10 cm thicknes

  G4double innerRadiusShield = 0*cm;
  G4double outerRadiusShield = 101*cm;
  G4double halfheightShield       = 5*cm;
  G4double startAngleShield  = 0*deg;
  G4double spanningAngleShield    = 360*deg;

  G4Tubs *aShield = new G4Tubs("aShield",
                               innerRadiusShield,
                               outerRadiusShield,
                               halfheightShield,
                               startAngleShield,
                               spanningAngleShield);

  G4Tubs *aShieldI1 = new G4Tubs("aShieldI1",
				 innerRadiusShield,
				 50*cm,
				 1*cm,
				 startAngleShield,
				 spanningAngleShield);

  
  // logical parallel cells


  G4LogicalVolume *aShield_logI1 = 
    new G4LogicalVolume(aShieldI1, fGalactic, "aShieldI1_log");

  fLogicalVolumeVector.push_back(aShield_logI1);  //ASO


  G4LogicalVolume *aShield_logI2 = 
    new G4LogicalVolume(aShieldI1, fGalactic, "aShieldI2_log");

  fLogicalVolumeVector.push_back(aShield_logI2);  //ASO

  // physical parallel cells

  G4String name;
  G4int i = 1;
  G4double startz = -85*cm; 
  for (i=1; i<=18; ++i) {
   
    name = fPVolumeStore.GetCellName(i);
    G4LogicalVolume *aShield_log = 
      new G4LogicalVolume(aShield, fGalactic, "aShield_log");

    fLogicalVolumeVector.push_back(aShield_log);  //ASO
    
    G4VPhysicalVolume *pvolIMinus = 
      new G4PVPlacement(0, 
			G4ThreeVector(0, 0, -0.5*halfheightShield),
			aShield_logI1, 
			name + "I1-", 
			aShield_log, 
			false, 
			i, true); //ASO
    //			0);
    
    G4VPhysicalVolume *pvolIPlus = 
      new G4PVPlacement(0, 
			G4ThreeVector(0, 0, +0.5*halfheightShield),
			aShield_logI2, 
			name + "I1+", 
			aShield_log, 
			false, 
			i, true); //ASO
    //			0);
    
    G4double pos_x = 0*cm;
    G4double pos_y = 0*cm;
    G4double pos_z = startz + (i-1) * (2*halfheightShield);
    G4VPhysicalVolume *pvol = 
      new G4PVPlacement(0, 
			G4ThreeVector(pos_x, pos_y, pos_z),
			aShield_log, 
			name, 
			worldCylinder_log, 
			false, 
			i, true); //ASO
    //			0);
    G4GeometryCell cell(*pvol, 0);
    fPVolumeStore.AddPVolume(cell);
    G4GeometryCell cellM(*pvolIMinus, 0);
    fPVolumeStore.AddPVolume(cellM);
    G4GeometryCell cellP(*pvolIPlus, 0);
    fPVolumeStore.AddPVolume(cellP);    
  }

  // filling the rest of the world volumr behind the concrete with
  // another slob which should get the same importance value as the 
  // last slob
  innerRadiusShield = 0*cm;
  //bug  outerRadiusShield = scaling*110*cm;
  outerRadiusShield = 101*cm;
  //bug  halfheightShield       = scaling*10*cm;
  halfheightShield       = 7.5*cm;
  startAngleShield  = 0*deg;
  spanningAngleShield    = 360*deg;

  G4Tubs *aRest = new G4Tubs("Rest",
			     innerRadiusShield,
			     outerRadiusShield,
			     halfheightShield,
			     startAngleShield,
			     spanningAngleShield);
  
  G4LogicalVolume *aRest_log = 
    new G4LogicalVolume(aRest, fGalactic, "aRest_log");

  fLogicalVolumeVector.push_back(aRest_log);  //ASO

  name = fPVolumeStore.GetCellName(19);
    
  G4double pos_x = 0*cm;
  G4double pos_y = 0*cm;
  //  G4double pos_z = 100*cm;
  G4double pos_z = 97.5*cm;
  G4VPhysicalVolume *pvol = 
    new G4PVPlacement(0, 
		      G4ThreeVector(pos_x, pos_y, pos_z),
		      aRest_log, 
		      name, 
		      worldCylinder_log, 
		      false, 
		      19, true); //ASO???
  //		      0);
  G4GeometryCell cell(*pvol, 0);
  fPVolumeStore.AddPVolume(cell);
  

  SetSensitive();


}

G4VPhysicalVolume &Tst33ParallelGeometry::GetWorldVolumeAddress() const{
  //  return *fWorldVolume;
  return *ghostWorld;
}

G4VPhysicalVolume *Tst33ParallelGeometry::GetWorldVolume() {
  //  return *fWorldVolume;
  return ghostWorld;
}


G4GeometryCell Tst33ParallelGeometry::
GetGeometryCell(G4int i, const G4String &nameExt) const {
  return fPVolumeStore.GetGeometryCell(i, nameExt);
}


void Tst33ParallelGeometry::SetSensitive(){

  G4cout << " CREATING Sensitive hit collections " << G4endl;
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


//   G4PSPopulation*   scorer2 = new G4PSPopulation(psName="Population");  
//   MFDet->RegisterPrimitive(scorer2);

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


// G4VUserParallelWorld &Tst33ParallelGeometry::GetParallelWorld() const{
//   return *fParallelWorld;
// }
