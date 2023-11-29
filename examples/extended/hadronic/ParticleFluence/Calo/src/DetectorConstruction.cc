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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4RotationMatrix.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "DetectorMessenger.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "PrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :  
  fVacuum( nullptr ), fIron( nullptr ), fCopper( nullptr ), fTungsten( nullptr ),
  fLead( nullptr ), fUranium( nullptr ), fPbWO4( nullptr ), fPolystyrene( nullptr ),
  fLiquidArgon( nullptr ), fSilicon( nullptr ), fQuartz( nullptr ), fBrass( nullptr ),
  fAluminium( nullptr ), fGraphite( nullptr ),
  fAbsorberMaterial( nullptr ), fActiveMaterial( nullptr ),
  fExperimentalHall_log( nullptr ), fExperimentalHall_phys( nullptr ),
  fLogicCalo( nullptr ), fPhysiCalo( nullptr ),
  fLogicModule( nullptr ), fPhysiModule( nullptr ),
  fLogicAbsorber( nullptr ), fPhysiAbsorber( nullptr ),
  fLogicActive( nullptr ), fPhysiActive( nullptr ),
  fFieldMgr( nullptr ), fUniformMagField( nullptr ), 
  fDetectorMessenger( nullptr ), 
  // Default values.  ***LOOKHERE***
  fIsCalHomogeneous( false ),    // Sampling calorimeter.
  fIsUnitInLambda( false ),      // Unit of length for the absorber total length.
  fAbsorberTotalLength( 2.0*CLHEP::m ), 
  fCalorimeterRadius( 1.0*CLHEP::m ), 
  fActiveLayerNumber( 50 ),
  fActiveLayerSize( 4.0*CLHEP::mm ),
  fIsRadiusUnitInLambda( false ), // Unit of length for the radius bin size.
  // Extra
  fCaloLength( 2.0*CLHEP::m ),
  // Scoring part
  fLogicScoringUpDown( nullptr ),
  fPhysiScoringUpstream( nullptr ),
  fPhysiScoringDownstream( nullptr ),
  fLogicScoringSide( nullptr ),
  fPhysiScoringSide( nullptr )
{
  fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  DefineMaterials();
  fAbsorberMaterial = fIron;
  fActiveMaterial = fPolystyrene;
  fDetectorMessenger = new DetectorMessenger( this );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete fDetectorMessenger;
  delete fUniformMagField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct() {
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials() { 
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, pressure, temperature, fractionmass;
  G4String name, symbol;
  G4int nel, natoms;

  //--- elements

  a = 1.01*g/mole;
  G4Element* elH = new G4Element( name="Hydrogen", symbol="H2", z=1., a );

  a = 2.01*g/mole;
  //G4Element* elD = new G4Element( name="Deuterium", symbol="D", z=1., a );

  a = 4.*g/mole;
  //G4Element* elHe = new G4Element( name="Helium", symbol="He", z=2., a );

  a = 6.94*g/mole;
  //G4Element* elLi = new G4Element( name="Lithium", symbol="Li", z=3., a );

  a = 9.01*g/mole;
  //G4Element* elBe = new G4Element( name="Berillium", symbol="Be", z=4., a );

  a = 12.01*g/mole;
  G4Element* elC = new G4Element( name="Carbon", symbol="C", z=6., a );

  a = 14.01*g/mole;
  G4Element* elN = new G4Element( name="Nitrogen", symbol="N2", z=7., a );

  a = 16.*g/mole;
  G4Element* elO = new G4Element( name="Oxygen", symbol="O2", z=8., a );

  a = 20.18*g/mole;
  //G4Element* elNe = new G4Element( name="Neon", symbol="Ne", z=10., a );

  a = 22.99*g/mole;
  //G4Element* elNa = new G4Element( name="Sodium", symbol="Na", z=11., a );

  a = 26.98*g/mole;
  //G4Element* elAl = new G4Element( name="Aluminium", symbol="Al", z=13., a );

  a = 28.085*g/mole;
  G4Element* elSi = new G4Element( name="Silicon", symbol="Si", z=14., a );

  a = 40.08*g/mole;
  //G4Element* elCa = new G4Element( name="Calcium", symbol="Ca", z=20., a );

  a = 55.850*g/mole;
  //G4Element* elFe = new G4Element( name="Iron", symbol="Fe", z=26., a );

  a = 63.54*g/mole;
  G4Element* elCu = new G4Element( name="Copper", symbol="Cu", z=29., a );

  a = 65.41*g/mole; 
  G4Element* elZn = new G4Element( name="Zinc", symbol="Zn", z=30., a );

  a = 183.85*g/mole;
  G4Element* elW = new G4Element( name="Tungstenm", symbol="W", z=74., a );

  a = 207.19*g/mole;
  G4Element* elPb = new G4Element( name="Lead", symbol="Pb", z=82., a );

  a = 238.03*g/mole;
  //G4Element* elU = new G4Element(name="Uranium", symbol="U", z=92., a);

  //--- simple materials

  // Iron has a  X0 = 1.7585 cm  and  lambda_I = 16.760 cm.   
  density = 7.87*g/cm3;
  a = 55.85*g/mole;
  fIron = new G4Material( name="Iron", z=26., a, density );

  // Copper has a  X0 = 1.4353 cm  and  lambda_I = 15.056 cm.   
  density = 8.96*g/cm3;
  a = 63.54*g/mole;
  fCopper = new G4Material( name="Copper", z=29., a, density );

  // Tungsten has a  X0 = 0.35 cm  and  lambda_I = 9.5855 cm. 
  density = 19.3*g/cm3;
  a = 183.85*g/mole;
  fTungsten = new G4Material( name="Tungsten", z=74., a, density );

  // Lead has a  X0 = 0.56120 cm  and  lambda_I = 17.092 cm.  
  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  fLead = new G4Material( name="Lead", z=82., a, density );

  // Uranium has a  X0 = 0.31662 cm  and  lambda_I = 10.501 cm.  
  density =  18.95*g/cm3;
  a = 238.03*g/mole;
  fUranium = new G4Material( name="Uranium", z=92., a, density );

  // Liquid Argon has a  X0 = 10.971 cm  and  lambda_I = 65.769 cm.  
  density = 1.4*g/cm3;
  a = 39.95*g/mole;
  fLiquidArgon = new G4Material( name="LiquidArgon", z=18., a, density );

  density = 0.002*g/cm3;
  a = 39.95*g/mole;
  //G4Material* ArgonGas = new G4Material( name="ArgonGas", z=18., a, density );

  // Silicon has a  X0 = 9.3688 cm  and  lambda_I = 46.5436 cm 
  density = 2.33*g/cm3;
  a = 28.085*g/mole;
  fSilicon = new G4Material( name="Silicon", z=14., a, density );
  
  // Aluminium has a  X0 = 8.8959 cm  and  lambda_I = 39.7184 cm
  density = 2.7*g/cm3;
  a = 26.98*g/mole;
  fAluminium = new G4Material( name="Aluminium", z=13., a, density );
  
  // Graphite has a  X0 = 19.3213 cm  and  lambda_I = 38.8235 cm
  density = 2.210*g/cm3;
  a = 12.0107*g/mole;
  fGraphite = new G4Material( name="Graphite", z=6., a, density );

  density = 8.96*g/cm3;
  a = 58.69*g/mole;
  //G4Material* Nickel = new G4Material( name="Nickel", z=28., a, density );

  //--- mixtures

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material( name="Air", density, nel=2 );
  Air->AddElement(elN, 0.7);
  Air->AddElement(elO, 0.3);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;  // From PhysicalConstants.h .
  fVacuum = new G4Material( name="Vacuum", density, nel=1,
                            kStateGas, temperature, pressure );
  fVacuum->AddMaterial( Air, fractionmass=1. );

  // Plastic scintillator tiles (used both in CMS hadron calorimeter
  // and ATLAS hadron barrel calorimeter): 
  //     X0 = 42.4 cm  and  lambda_I = 79.360 cm.  
  density = 1.032*g/cm3;
  fPolystyrene = new G4Material( name="Polystyrene", density, nel=2 );
  fPolystyrene->AddElement( elC, natoms=19 );
  fPolystyrene->AddElement( elH, natoms=21 );

  // PbWO4 CMS crystals. It has a  X0 = 0.89 cm  and  lambda_I = 22.4 cm. 
  density = 8.28*g/cm3;
  fPbWO4 = new G4Material( name="PbWO4", density, nel=3 );
  fPbWO4->AddElement( elPb, natoms=1 );
  fPbWO4->AddElement( elW,  natoms=1 );
  fPbWO4->AddElement( elO,  natoms=4 );

  fQuartz = new G4Material( name="Quartz", density=2.200*g/cm3, nel=2 );
  fQuartz->AddElement( elSi, 1 );
  fQuartz->AddElement( elO , 2 );

  fBrass = new G4Material( name="Brass", density=8.6*g/cm3, nel=2 );
  fBrass->AddElement( elCu, 0.7 );
  fBrass->AddElement( elZn, 0.3 );  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter() {
  if ( ! AreParametersOK() ) {
    G4cout << " DetectorConstruction::ConstructCalorimeter() : ***ERROR*** "
           << G4endl << "\t PARAMETERS NOT WELL-DEFINED! GEOMETRY UNCHANGED."
           << G4endl;
    return fExperimentalHall_phys;
  }

  // Clean old geometry, if any.
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4double lambda = 0.0; // G4double X0 = 0.0;
  if ( fIsUnitInLambda ) {
    if ( fAbsorberMaterial == fIron ) {
      lambda = 16.760*cm; // X0     = 1.7585*cm; 
    } else if ( fAbsorberMaterial == fCopper ) {
      lambda = 15.056*cm; // X0     = 1.4353*cm; 
    } else if ( fAbsorberMaterial == fBrass ) {  
      lambda = 15.056*cm; // Lack of PDG data: I am assuming the same as Copper.  // X0=1.4353*cm
    } else if ( fAbsorberMaterial == fTungsten ) {
      lambda = 9.5855*cm; // X0     = 0.35*cm; 
    } else if ( fAbsorberMaterial == fLead ) {
      lambda = 17.092*cm; // X0     = 0.56120*cm; 
    } else if ( fAbsorberMaterial == fPbWO4 ) {
      lambda = 22.4*cm;   // X0     = 0.89*cm; 
    } else if ( fAbsorberMaterial == fUranium ) {
      lambda = 10.501*cm; // X0     = 0.31662*cm; 
    } else if ( fAbsorberMaterial == fGraphite ) {
      lambda = 38.82*cm;  // X0     = 19.32*cm; 
    } else {
      std::cout << "ERROR: absorber material not recognized" << std::endl;
    }
  }

  //------------------- volumes --------------------------

  G4double absorberTotalLength = fAbsorberTotalLength;
  G4double calorimeterRadius = fCalorimeterRadius;
  if ( fIsUnitInLambda ) {
    absorberTotalLength *= lambda; 
    calorimeterRadius *= lambda; 
  }

  // --- experimental hall (world volume)    ***LOOKHERE***
  //     beam line along the Z-axis
  G4double expHall_x = 10.0*m;  // half dimension along x 
  G4double expHall_y = 10.0*m;  // half dimension along y
  G4double expHall_z = 10.0*m;  // half dimension along z

  G4Box* experimentalHall_box = new G4Box( "expHall_box", expHall_x, expHall_y, expHall_z );
  fExperimentalHall_log = new G4LogicalVolume( experimentalHall_box,  // solid 
                                               fVacuum,               // material
                                               "expHall_log",         // name
                                               0,                     // field manager
                                               0,                     // sensitive detector
                                               0 );                   // user limits
  fExperimentalHall_phys = new G4PVPlacement( 0,                      // rotation
                                              G4ThreeVector(),        // translation
                                              "expHall",              // name
                                              fExperimentalHall_log,  // logical volume
                                              0,                      // mother physical volume
                                              false,                  // boolean operation
                                              0 );                    // copy number
  
  // --- Detector
  // The idea is to use Replica placement. 
  // To do that, we have to define two extra volumes: the "calorimeter" volume 
  // and the "module". The former, which has the world as its mother volume, 
  // is the mother of the module volume. The calorimeter volume is completely 
  // filled by a number (theActiveLayerNumber) of replicas of the module volume. 
  // A module volume, in its turn, is the mother volume of the absorber layer + 
  // active layer. 

  //            --- absorber layer : logical
  G4double zAbsorber = absorberTotalLength / static_cast< double >( fActiveLayerNumber );
  // In the case of homogenous calorimeter the "active" part must be
  // subtracted because it is made of the same material
  // (the material of the "active" part is set to be the same as
  //  the aborber).
  if ( fIsCalHomogeneous ) { 
    fActiveMaterial = fAbsorberMaterial;
    zAbsorber -= fActiveLayerSize; 
  }
  zAbsorber /= 2.0;  // half dimension along z
  G4Tubs* solidAbsorber = new G4Tubs( "solidAbsorber",      // name
                                      0.0,                  // inner radius
                                      calorimeterRadius,    // outer radius
                                      zAbsorber,            // half cylinder length in z
                                      0.0,                  // starting phi angle in rad
                                      2.0*pi );             // final phi angle in rad
  fLogicAbsorber = new G4LogicalVolume( solidAbsorber,      // solid 
                                        fAbsorberMaterial,  // material
                                        "logicAbsorber",    // name
                                        0,                  // field manager
                                        0,                  // sensitive detector
                                        0 );                // user limits

  //            --- active layer : logical
  G4double zActive = fActiveLayerSize / 2.0;  // half dimension along z 
  G4Tubs* solidActive = new G4Tubs( "solidActive",      // name
                                    0.0,                // inner radius
                                    calorimeterRadius,  // outer radius
                                    zActive,            // half cylinder length in z
                                    0.0,                // starting phi angle in rad
                                    2.0*pi );           // final phi angle in rad
  fLogicActive = new G4LogicalVolume( solidActive,      // solid 
                                      fActiveMaterial,  // material
                                      "logicActive",    // name
                                      0,                // field manager
                                      0,                // sensitive detector
                                      0 );              // user limits

  //        --- module : logical
  G4double zModule = zAbsorber + zActive;  // half dimension along z 
  G4Tubs* solidModule = new G4Tubs( "solidModule",      // name
                                    0.0,                // inner radius
                                    calorimeterRadius,  // outer radius
                                    zModule,            // half cylinder length in z
                                    0.0,                // starting phi angle in rad
                                    2.0*pi );           // final phi angle in rad
  fLogicModule = new G4LogicalVolume( solidModule,      // solid 
                                      fLead,            // material, it does NOT matter
                                      "logicModule",    // name
                                      0,                // field manager
                                      0,                // sensitive detector
                                      0 );              // user limits

  //    --- calorimeter : logical
  G4int numberOfModules = fActiveLayerNumber;
  G4double zCalo = numberOfModules*zModule;  // half dimension along z
  fCaloLength = 2.0*zCalo;
  G4Tubs* solidCalo = new G4Tubs( "solidCalo",        // name
                                  0.0,                // inner radius
                                  calorimeterRadius,  // outer radius
                                  zCalo,              // half cylinder length in z
                                  0.0,                // starting phi angle in rad
                                  2.0*pi );           // final phi angle in rad
  fLogicCalo = new G4LogicalVolume( solidCalo,        // solid 
                                    fLead,            // material, it does NOT matter
                                    "logicCalo",      // name
                                    0,                // field manager
                                    0,                // sensitive detector
                                    0 );              // user limits

  //            --- absorber layer : physical
  G4double zpos = - zActive;
  fPhysiAbsorber = new G4PVPlacement( 0,                        // rotation
                                      G4ThreeVector(0,0,zpos),  // translation
                                      fLogicAbsorber,           // logical volume
                                      "physiAbsorber",          // name
                                      fLogicModule,             // mother logical volume
                                      false,                    // boolean operation
                                      1000 );                   // copy number

  //            --- active layer : physical
  zpos += zAbsorber + zActive;
  fPhysiActive = new G4PVPlacement( 0,                          // rotation
                                    G4ThreeVector(0,0,zpos),    // translation
                                    fLogicActive,               // logical volume
                                    "physiActive",              // name
                                    fLogicModule,               // mother logical volume
                                    false,                      // boolean operation
                                    2000 );                     // copy number

  //        --- module : physical (using replica)
  fPhysiModule = new G4PVReplica( "Calo",                   // name
                                  fLogicModule,             // logical volume
                                  fLogicCalo,               // mother logical volume
                                  kZAxis,                   // axis of replication
                                  numberOfModules,          // number of replica
                                  2*(zAbsorber+zActive) );  // (full) width of replica

  //    --- calorimeter : physical
  fPhysiCalo = new G4PVPlacement( 0,                        // rotation
                                  G4ThreeVector(),          // translation
                                  "physiCalo",              // its name
                                  fLogicCalo,               // logical volume
                                  fExperimentalHall_phys,   // mother physical volume
                                  false,                    // boolean operation
                                  100 );                    // copy number

  // Three scoring volumes: one thin layer downstream of the calorimeter ("down")
  //                        one thin layer surrounding (lateral) of the calorimeter ("side")
  //                        one thin layer upstream of the calorimeter ("up")
  G4Tubs* solidScoringUpDown = new G4Tubs( "solidScoringUpDown",    // name
                                           0.0,                     // inner radius
                                           calorimeterRadius,       // outer radius
                                           0.5*fScoringThickness,   // half cylinder length in z
                                           0.0,                     // starting phi angle in rad
                                           2.0*pi );                // final phi angle in rad
  fLogicScoringUpDown = new G4LogicalVolume( solidScoringUpDown,    // solid
                                             fVacuum,               // material
                                             "logicScoringUpDown",  // name
                                             0,                     // field manager
                                             0,                     // sensitive detector
                                             0 );                   // user limits
  G4double zScoringUpDown = 0.5*(fCaloLength + fScoringThickness);
  fPhysiScoringUpstream = new G4PVPlacement( 0,                          // rotation
                                             G4ThreeVector( 0.0, 0.0, -zScoringUpDown ),
                                                                         // translation
                                             "physiScoringUpstream",     // name
                                             fLogicScoringUpDown,        // logical volume
                                             fExperimentalHall_phys,     // mother physical volume
                                             false,                      // boolean operation
                                             0 );                        // copy number
  fPhysiScoringDownstream = new G4PVPlacement( 0,                        // rotation
                                               G4ThreeVector( 0.0, 0.0, zScoringUpDown ),
                                                                         // translation
                                               "physiScoringDownstream", // name
                                               fLogicScoringUpDown,      // logical volume
                                               fExperimentalHall_phys,   // mother physical volume
                                               false,                    // boolean operation
                                               0 );                      // copy number
  
  G4Tubs* solidScoringSide = new G4Tubs( "solidScoringSide",    // name
                                         calorimeterRadius,     // inner radius
                                         calorimeterRadius + fScoringThickness,  // outer radius
                                         0.5*fCaloLength,       // half cylinder length in z
                                         0.0,                   // starting phi angle in rad
                                         2.0*pi );              // final phi angle in rad
  fLogicScoringSide = new G4LogicalVolume( solidScoringSide,    // solid
                                           fVacuum,             // material
                                           "logicScoringSide",  // name
                                           0,                   // field manager
                                           0,                   // sensitive detector
                                           0 );                 // user limits
  fPhysiScoringSide = new G4PVPlacement( 0,                       // rotation
                                         G4ThreeVector( 0.0, 0.0, 0.0 ),  // translation
                                         "physiScoringSide",      // name
                                         fLogicScoringSide,       // logical volume
                                         fExperimentalHall_phys,  // mother physical volume
                                         false,                   // boolean operation
                                         0 );                     // copy number
  
  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DetectorConstruction::AreParametersOK() {
  bool isOk = true;
  if ( ! fAbsorberMaterial ) {
    isOk = false;
    G4cout << " DetectorConstruction::AreParametersOK() : UNDEFINED absorber material" << G4endl;
  }
  if ( ! fActiveMaterial ) {
    isOk = false;
    G4cout << " DetectorConstruction::AreParametersOK() : UNDEFINED active material" << G4endl;
  }
  if ( fAbsorberTotalLength <= 0.0 ) {
    isOk = false;
    G4cout << " DetectorConstruction::AreParametersOK() : fAbsorberTotalLength = "
           << fAbsorberTotalLength << G4endl;
  }
  if ( fCalorimeterRadius <= 0.0 ) {
    isOk = false;
    G4cout << " DetectorConstruction::AreParametersOK() : fCalorimeterRadius = "
           << fCalorimeterRadius << G4endl;
  }
  if ( fActiveLayerNumber <= 0 ) {
    isOk = false;
    G4cout << " DetectorConstruction::AreParametersOK() : fActiveLayerNumber = "
           << fActiveLayerNumber << G4endl;
  }
  if ( fActiveLayerSize <= 0.0 ) {
    isOk = false;
    G4cout << " DetectorConstruction::AreParametersOK() : fActiveLayerSize = "
           << fActiveLayerSize << G4endl;
  }
  return isOk;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField( const G4double fieldValue ) {
  if ( fUniformMagField ) {
    delete fUniformMagField;
  }
  if ( std::abs( fieldValue ) > 0.0 ) {
    // Apply a global uniform magnetic field along the Y axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.
    fUniformMagField = new G4UniformMagField( G4ThreeVector( 0.0, fieldValue, 0.0 ) );
    fFieldMgr->SetDetectorField( fUniformMagField );
    fFieldMgr->CreateChordFinder( fUniformMagField );
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetAbsorberMaterial( const G4String name ) {
  if        ( name == "Fe" || name == "Iron" || name == "iron" ) { 
    fAbsorberMaterial = fIron;
  } else if ( name == "Cu" || name == "Copper" || name == "copper" ) { 
    fAbsorberMaterial = fCopper;
  } else if ( name == "Brass" || name == "brass" ) { 
    fAbsorberMaterial = fBrass;
  } else if ( name == "Pb" || name == "Lead" || name == "lead" ) { 
    fAbsorberMaterial = fLead;
  } else if ( name == "PbWO4" ) {
    fAbsorberMaterial = fPbWO4;
  } else if ( name == "W" || name == "Tungsten" || name == "tungsten" ) { 
    fAbsorberMaterial = fTungsten;
  } else if ( name == "U" || name == "Uranium" || name == "uranium" ) { 
    fAbsorberMaterial = fUranium;
  } else if ( name == "C" || name == "Graphite" || name == "graphite" ) { 
    fAbsorberMaterial = fGraphite;
  } else {
    G4cout << G4endl << G4endl
        << "WARNING: the name of the material has not been recognized!" << G4endl
        << "     ===> the default  * Iron *  will be used." << G4endl << G4endl;  
    fAbsorberMaterial = fIron;
  }  
  fLogicAbsorber->SetMaterial( fAbsorberMaterial );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetActiveMaterial( const G4String name ) {
  if        ( name == "Scintillator" || name == "scintillator" ) { 
    fActiveMaterial = fPolystyrene;
  } else if ( name == "LAr" || name == "LiquidArgon" || name == "liquidArgon" ) { 
    fActiveMaterial = fLiquidArgon;
  } else if ( name == "PbWO4" ) {
    fActiveMaterial = fPbWO4;
  } else if ( name == "Si" || name == "Silicon" || name == "silicon" ) { 
    fActiveMaterial = fSilicon;
  } else if ( name == "Quartz" || name == "quartz" ) { 
    fActiveMaterial = fQuartz;
  } else if ( name == "C" || name == "Graphite" || name == "graphite" ) { 
    fActiveMaterial = fGraphite;
  } else {
    G4cout << G4endl << G4endl
        << "WARNING: the name of the material has not been recognized!" << G4endl
        << "     ===> the default  * Scintillator *  will be used." << G4endl << G4endl;
    fActiveMaterial = fPolystyrene;
  }  
  fLogicActive->SetMaterial( fActiveMaterial );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::UpdateGeometry() {
  //G4RunManager::GetRunManager()->DefineWorldVolume( ConstructCalorimeter() );
  G4RunManager::GetRunManager()->ReinitializeGeometry();
  PrintParameters();
  // Update also the position of the gun
  const PrimaryGeneratorAction* pPrimaryAction = 
    dynamic_cast< const PrimaryGeneratorAction* >(
      G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction() );
  if ( pPrimaryAction ) pPrimaryAction->SetGunPosition();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() {
  G4cout << G4endl << G4endl << " ------  DetectorConstruction::PrintParameters() ------ "
         << G4endl
         << " Absorber Material = ";
  if ( fAbsorberMaterial ) {
    G4cout << fAbsorberMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << " Active Material   = ";
  if ( fActiveMaterial ) {
    G4cout << fActiveMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << " Is the Calorimeter Homogeneous ? " << fIsCalHomogeneous;
  G4cout << G4endl << " Is the Unit in Lambda ? " << fIsUnitInLambda;
  G4cout << G4endl << " Absorber Total Length = ";
  if ( fIsUnitInLambda ) {
    G4cout << fAbsorberTotalLength << "  lambdas";
  } else {
    G4cout << fAbsorberTotalLength / m << " m";
  }
  G4cout << G4endl << " Calorimeter Radius = ";
  if ( fIsUnitInLambda ) {
    G4cout << fCalorimeterRadius << "  lambdas";
  } else {
    G4cout << fCalorimeterRadius / m << " m";
  }
  G4cout << G4endl << " Active Layer Number   = " << fActiveLayerNumber;
  G4cout << G4endl << " Active Layer Size     = " << fActiveLayerSize/mm << " mm";
  G4cout << G4endl << " Is the Radius Unit in Lambda ? " << fIsRadiusUnitInLambda;
  G4cout << G4endl << " Radius Bin Size       = ";
  G4cout << G4endl << " Magnetic field [T]    = ";
  if ( fUniformMagField ) {
    G4cout << fUniformMagField->GetConstantFieldValue() / tesla;
  } else {
    G4cout << "(0,0,0)";
  }

  G4cout << G4endl << " -------------------------------------------------------- " << G4endl
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
