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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISADetectorConstruction class                                   *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************



#include "LISADetectorConstruction.hh"
#include "LISADetectorMaterials.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISADetectorConstruction::LISADetectorConstruction() {;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
LISADetectorConstruction::~LISADetectorConstruction() {;}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* LISADetectorConstruction::Construct() {


  // material definitions
  ConstructMaterials();


  // build it
  return ConstructDetector();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* LISADetectorConstruction::ConstructDetector() {


  // Edge and solid colour attributes
#include "LISAColours.icc"



  //***************************************************************************
  // world
  //***************************************************************************

  G4double wld_len = 3000.*mm;

  G4Tubs* truewld_sol  =
    new G4Tubs("truewld_box", 0., 0.5*wld_len, 0.5*wld_len, 0., 360.*deg);
  G4LogicalVolume* truewld_log  = 
    new G4LogicalVolume(truewld_sol, vacuum, "truewld_log");
  G4VPhysicalVolume* truewld_phys = new G4PVPlacement(0, G4ThreeVector(), 
    "truewld_phys", truewld_log, NULL, false, 0);
  truewld_log->SetVisAttributes(G4VisAttributes::Invisible);

  // allow spacecraft rotation
  G4RotationMatrix wld_rot; wld_rot.rotateZ(0.*deg);
  G4Tubs* wld_sol  =
    new G4Tubs("wld_box", 0., 0.49*wld_len, 0.5*wld_len, 0., 360.*deg);
  G4LogicalVolume* wld_log  = 
    new G4LogicalVolume(wld_sol, vacuum, "wld_log");
  G4VPhysicalVolume* wld_phys = new G4PVPlacement(G4Transform3D(wld_rot, 
    G4ThreeVector()), "wld_phys", wld_log, truewld_phys, false, 0);
  //  wld_log->SetVisAttributes(white_vat);
  wld_log->SetVisAttributes(G4VisAttributes::Invisible);


  // Probe volume
  // uncomment code in LISASteppingAction.cc to count hits
  //   G4Sphere* probe_sol = 
  //     new G4Sphere("probe_sol", 0., 100.*mm, 0., 360.*deg, 0., 180.*deg);
  //   G4LogicalVolume* probe_log  = 
  //     new G4LogicalVolume(probe_sol, vacuum, "probe_log");
  //   G4VPhysicalVolume* probe_phys = 
  //     new G4PVPlacement(0, G4ThreeVector(0*mm,0*mm,0*mm),
  // 	 "probe_phys", probe_log, wld_phys, false, 0);
  //   probe_log->SetVisAttributes(sol_white_vat);



  //**************************************************************************
  // Science Module Structure (SMS)
  //**************************************************************************
  //    Primary Structure
  //    Lower Deck
  //    Upper Deck
  //    Thermal Shield
  //    Solar Array 
  //    Optical Surface Reflectors (OSR)
  //    Radiator Panels
  //**************************************************************************

  // position of science module spacecraft in world volume
  G4ThreeVector spacecraft_pos(0.*mm,0.*mm,0.*mm);

#include "LISAScienceModuleStructures.icc"



  //***************************************************************************
  // Interferometer Assembly
  //***************************************************************************
  //   Payload Shield (Y Tube)
  //   Telescope Light Shields
  //   Telescope Electronics Mounting
  //   Telescope Mirrors Mounting
  //   Telescope Actuator Mechanism
  //   Optical Bench Mounting
  //   Optical Bench
  //***************************************************************************


  // position of payload shields in spacecraft
  G4double YTube_xoff = -800.*mm;

  // position of interferometer in Y-tubes
  G4RotationMatrix tel_rot1; tel_rot1.rotateY(210.*deg);
  G4ThreeVector tel_pos1(+319.*mm, 0.*mm, 920.*mm);
  G4RotationMatrix tel_rot2; tel_rot2.rotateY(150.*deg);
  G4ThreeVector tel_pos2(-319.*mm, 0.*mm, 920.*mm);

  // position of optical bench along Y-tube 
  G4double OpticalBench_off = 0.0*mm;

#include "LISAInterferometerAssembly.icc"



  //***************************************************************************
  // Sensor Vacuum Housing
  // Modified to contain LTP Inertial Sensor and Caging Mechanism
  //***************************************************************************

  G4RotationMatrix IS_rot; IS_rot.rotateX(90.*deg); IS_rot.rotateY(90.*deg);

#include "LISASensorHousing.icc"



  //***************************************************************************
  // LTP Caging Mechanism
  //***************************************************************************

  //#include "LISACagingMechanism.icc"




  //***************************************************************************
  // Inertial Sensors:  YZ-Injection, 46 mm Test Mass
  // Design adopted for LTP/SMART-2: Report LTP-RT-CGS-001, issue 2
  //***************************************************************************

#include "LISAInertialSensor.icc"

  // Sensor Region (cuts 250 eV)
  G4Region* ISensor = new G4Region(G4String("sensor"));
  cage_o_log->SetRegion(ISensor);
  ISensor->AddRootLogicalVolume(cage_o_log);




  //***************************************************************************
  // Electronics Boxes
  //***************************************************************************

#include "LISAElectronicsBoxes.icc"




  //***************************************************************************
  // Support Systems
  //***************************************************************************
  //   Star Trackers
  //   FEEP Thrusters
  //   Communications Antennas
  //***************************************************************************

#include "LISASupportSystems.icc"



  // ......................................................................
  // attach user limits ...................................................

  // reduce step size in electrodes
  // goldplating_log->SetUserLimits (new G4UserLimits(30.*nanometer));

  // return
  return truewld_phys;

}
