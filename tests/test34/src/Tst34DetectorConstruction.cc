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
#include <iostream>

// Geant4 Classes
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4Colour.hh"

// User Classes
#include "Tst34DetectorConstruction.hh"
#include "Tst34SensitiveDetector.hh"
#include "Tst34MaterialManager.hh"

// Fast simulation
#include "GFlashHomoShowerParameterisation.hh"
#include "G4FastSimulationManager.hh"
#include "GFlashShowerModel.hh"
#include "GFlashHitMaker.hh"
#include "GFlashParticleBounds.hh"

Tst34DetectorConstruction::Tst34DetectorConstruction()
 : m_experimentalHall_log(0), 
   m_calo_log(0), 
   m_experimentalHall_phys(0), 
   m_calo_phys(0)
{
  G4cout << "Tst34DetectorConstruction::Detector constructor" << G4endl;    
  
  // Simplified `CMS-like` PbWO4 crystal calorimeter  
  m_calo_xside=31*cm;
  m_calo_yside=31*cm;
  m_calo_zside=24*cm; 
  
  // GlashStuff
  m_theParticleBounds  = new GFlashParticleBounds(); // Energy Cuts to kill particles
  m_theHMaker          = new GFlashHitMaker();       // Makes the EnergieSpots
}

Tst34DetectorConstruction::~Tst34DetectorConstruction()
{ 
  delete m_theParametrisation;
  delete m_theParticleBounds;
  delete m_theFastShowerModel;
}

G4VPhysicalVolume* Tst34DetectorConstruction::Construct()
{
  // -------- Definitions of Solids, Logical Volumes, Physical Volumes --------
  G4String mat= "PbWO4";     
  G4cout << "Defining the materials" << G4endl;
  Tst34MaterialManager *matManager=Tst34MaterialManager::GetMaterialManager();

  /*******************************
  * The Experimental Hall       *
  *******************************/
  m_experimentalHall_x=1000.*cm;
  m_experimentalHall_y=1000.*cm;
  m_experimentalHall_z=1000.*cm;

  m_experimentalHall_box = new G4Box("expHall_box",         // World Volume
                                     m_experimentalHall_x,  // x size
                                     m_experimentalHall_y,  // y size 
                                     m_experimentalHall_z); // z size

  m_experimentalHall_log = new G4LogicalVolume(m_experimentalHall_box,
                                               matManager->getMaterial("Air"),
                                               "expHall_log",   // its name 
                                               0,     //opt: fieldManager
                                               0,     //opt: SensitiveDetector 
                                               0);    //opt: UserLimits                       
  m_experimentalHall_phys = new G4PVPlacement(0,  // no rotation       
                                G4ThreeVector(),  // at (0,0,0)
                                "expHall",        // its name 
                                m_experimentalHall_log, // its logical volume
                                0,                      // its mother  volume
                                false,                  // no boolean operation
                                0);                     // copy number

  //------------------------------ 
  // Calorimeter segments
  //------------------------------
  // Simplified `CMS-like` PbWO4 crystal calorimeter  

  G4Box *calo_box= new G4Box("CMS calorimeter",     // its name
                             m_calo_xside/2.,       // size   
                             m_calo_yside/2.,
                             m_calo_zside/2.);
  m_calo_log = new G4LogicalVolume(calo_box,
                                   matManager->getMaterial("PbWO4"),
                                   "calo log",  // its name
                                   0,           // opt: fieldManager
                                   0,           // opt: SensitiveDetector 
                                   0);          // opt: UserLimit
  G4double Xpos = 0.0;
  G4double Ypos = 0.0;
  G4double Zpos = 100.0*cm;

  m_calo_phys = new G4PVPlacement(0,                             // no rotation
                                  G4ThreeVector(Xpos,Ypos,Zpos), // at (0,0,0)
                                  m_calo_log,             // its logical volume     
                                  "calorimeter",          // its name
                                  m_experimentalHall_log, // its mother volume
                                  false,                  // no boolean
                                  1);                     // visibility

  // Sensitive Detector part
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  Tst34SensitiveDetector* CaloSD=
  new Tst34SensitiveDetector("Calorimeter",this);
  SDman->AddNewDetector(CaloSD);
 
  m_experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* CaloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));

  m_calo_log->SetVisAttributes(CaloVisAtt);
  m_calo_log->SetSensitiveDetector(CaloSD);

  /**********************************************
  * Initializing shower modell
  ***********************************************/
  G4cout << "Shower parameterization" << G4endl;
  m_theFastShowerModel = new GFlashShowerModel("fastShowerModel",m_calo_log);
  m_theParametrisation =
    new GFlashHomoShowerParameterisation(matManager->getMaterial(mat));
  m_theFastShowerModel->SetParametrisation(*m_theParametrisation);
  m_theFastShowerModel->SetParticleBounds(*m_theParticleBounds) ;
  m_theFastShowerModel->SetHitMaker(*m_theHMaker); 
  G4cout << "end shower parameterization" << G4endl;
  /**********************************************/

  return m_experimentalHall_phys;
}
