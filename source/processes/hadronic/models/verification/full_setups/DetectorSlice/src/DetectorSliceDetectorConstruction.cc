#include "DetectorSliceDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"

#include "G4PVReplica.hh"

#include "G4TransportationManager.hh"
      
#include "DetectorSliceDetectorMessenger.hh"

#include "DetectorSliceSensitiveEmCalo.hh"
#include "DetectorSliceSensitiveHadCalo.hh"
#include "G4SDManager.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
           

DetectorSliceDetectorConstruction::DetectorSliceDetectorConstruction() :  
  Vacuum( 0 ), Iron( 0 ), Copper( 0 ), Tungsten( 0 ), Lead( 0 ), Uranium( 0 ), 
  PbWO4( 0 ), Polystyrene( 0 ), LiquidArgon( 0 ), Silicon( 0 ), Quartz( 0 ),
  theTrackerMaterial(0), theEmAbsorberMaterial( 0 ), theEmActiveMaterial( 0 ),
  theHadAbsorberMaterial( 0 ), theHadActiveMaterial( 0 ), theMuonMaterial( 0 ),
  experimentalHall_log( 0 ), experimentalHall_phys( 0 ),
  logicTracker( 0 ), physiTracker( 0 ),
  logicEmCalo( 0 ), physiEmCalo( 0 ),
  logicEmModule( 0 ), physiEmModule( 0 ),
  logicEmAbsorber( 0 ), physiEmAbsorber( 0 ),
  logicEmActive( 0 ), physiEmActive( 0 ),
  logicHadCalo( 0 ), physiHadCalo( 0 ),
  logicHadModule( 0 ), physiHadModule( 0 ),
  logicHadAbsorber( 0 ), physiHadAbsorber( 0 ),
  logicHadActive( 0 ), physiHadActive( 0 ),
  logicMuon( 0 ), physiMuon( 0 ),
  detectorMessenger( 0 ), 
  theSensitiveEmCalorimeter( 0 ), theSensitiveHadCalorimeter( 0 ), 
  theVisAttAbsorber( 0 ), theVisAttActive( 0 ),

  // Default values.  ***LOOKHERE***
  theIsEmCalHomogeneous( false ),    // Sampling EM calorimeter.
  theIsHadCalHomogeneous( false ),   // Sampling HAD calorimeter.
  theTrackerLength( 1.0*cm ), 
  theEmAbsorberTotalLength( 20.0*cm ), 
  theHadAbsorberTotalLength( 2.0*m ), 
  theMuonLength( 1.0*m ), 
  theDetectorRadius( 1.0*m ), 
  theEmActiveLayerNumber( 50 ),
  theHadActiveLayerNumber( 50 ),
  theEmActiveLayerSize( 1.0*mm ),
  theHadActiveLayerSize( 4.0*mm )
{
  //G4cout << " BEGIN  DetectorSliceDetectorConstruction::DetectorSliceDetectorConstruction()" << G4endl; //***DEBUG***

  // materials
  DefineMaterials();

  theTrackerMaterial     = Silicon;
  theEmAbsorberMaterial  = Lead;
  theEmActiveMaterial    = LiquidArgon;
  theHadAbsorberMaterial = Iron;
  theHadActiveMaterial   = Polystyrene;
  theMuonMaterial        = Iron;

  detectorMessenger = new DetectorSliceDetectorMessenger( this );

  //G4cout << " END  DetectorSliceDetectorConstruction::DetectorSliceDetectorConstruction()" << G4endl; //***DEBUG***
}


DetectorSliceDetectorConstruction::~DetectorSliceDetectorConstruction() {
  delete detectorMessenger;
}


G4VPhysicalVolume* DetectorSliceDetectorConstruction::Construct() {
  //G4cout << " BEGIN  DetectorSliceDetectorConstruction::Construct()" << G4endl; //***DEBUG***

  return ConstructCalorimeter();
}


void DetectorSliceDetectorConstruction::DefineMaterials() { 

  //G4cout << " BEGIN  DetectorSliceDetectorConstruction::DefineMaterials()" 
  //       << G4endl; //***DEBUG***

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
  //G4Element* elCu = new G4Element( name="Copper", symbol="Cu", z=29., a );

  a = 183.85*g/mole;
  G4Element* elW = new G4Element( name="Tungstenm", symbol="W", z=74., a );

  a = 207.19*g/mole;
  G4Element* elPb = new G4Element( name="Lead", symbol="Pb", z=82., a );

  a = 238.03*g/mole;
  //G4Element* elU = new G4Element(name="Uranium", symbol="U", z=92., a);

  //--- simple materials

  density = 2.7*g/cm3;
  a = 26.98*g/mole;
  //G4Material* Aluminium = new G4Material( name="Aluminium", z=13., a, density );
  
  // Iron has a  X0 = 1.7585 cm  and  lambda_I = 16.760 cm.   
  density = 7.87*g/cm3;
  a = 55.85*g/mole;
  Iron = new G4Material( name="Iron", z=26., a, density );

  // Copper has a  X0 = 1.4353 cm  and  lambda_I = 15.056 cm.   
  density = 8.96*g/cm3;
  a = 63.54*g/mole;
  Copper = new G4Material( name="Copper", z=29., a, density );

  // Tungsten has a  X0 = 0.35 cm  and  lambda_I = 9.5855 cm. 
  density = 19.3*g/cm3;
  a = 183.85*g/mole;
  Tungsten = new G4Material( name="Tungsten", z=74., a, density );

  // Lead has a  X0 = 0.56120 cm  and  lambda_I = 17.092 cm.  
  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  Lead = new G4Material( name="Lead", z=82., a, density );

  // Uranium has a  X0 = 0.31662 cm  and  lambda_I = 10.501 cm.  
  density =  18.95*g/cm3;
  a = 238.03*g/mole;
  Uranium = new G4Material( name="Uranium", z=92., a, density );

  // Liquid Argon has a  X0 = 10.971 cm  and  lambda_I = 65.769 cm.  
  density = 1.4*g/cm3;
  a = 39.95*g/mole;
  LiquidArgon = new G4Material( name="LiquidArgon", z=18., a, density );

  density = 0.002*g/cm3;
  a = 39.95*g/mole;
  //G4Material* ArgonGas = new G4Material( name="ArgonGas", z=18., a, density );

  density = 2.33*g/cm3;
  a = 28.085*g/mole;
  Silicon = new G4Material( name="Silicon", z=14., a, density );
  
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
  Vacuum = new G4Material( name="Vacuum", density, nel=1,
			   kStateGas, temperature, pressure );
  Vacuum->AddMaterial( Air, fractionmass=1. );

  // Plastic scintillator tiles (used both in CMS hadron calorimeter
  // and ATLAS hadron barrel calorimeter): 
  //     X0 = 42.4 cm  and  lambda_I = 79.360 cm.  
  density = 1.032*g/cm3;
  Polystyrene = new G4Material( name="Polystyrene", density, nel=2 );
  Polystyrene->AddElement( elC, natoms=19 );
  Polystyrene->AddElement( elH, natoms=21 );

  // PbWO4 CMS crystals. It has a  X0 = 0.89 cm  and  lambda_I = 22.4 cm. 
  density = 8.28*g/cm3;
  PbWO4 = new G4Material( name="PbWO4", density, nel=3 );
  PbWO4->AddElement( elPb, natoms=1 );
  PbWO4->AddElement( elW,  natoms=1 );
  PbWO4->AddElement( elO,  natoms=4 );

  Quartz = new G4Material( name="Quartz", density=2.200*g/cm3, nel=2 );
  Quartz->AddElement( elSi, 1 );
  Quartz->AddElement( elO , 2 );

  //G4cout << " END  DetectorSliceDetectorConstruction::DefineMaterials()" 
  //       << G4endl; //***DEBUG***

}


G4VPhysicalVolume* DetectorSliceDetectorConstruction::ConstructCalorimeter() {

  //G4cout << " BEGIN  DetectorSliceDetectorConstruction::ConstructCalorimeter()" 
  //       << G4endl; //***DEBUG***

  if ( ! areParametersOK() ) {
    G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : ***ERROR*** "
           << G4endl << "\t PARAMETERS NOT WELL-DEFINED! GEOMETRY UNCHANGED."
           << G4endl;
    return experimentalHall_phys;
  }

  // Clean old geometry, if any.
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //------------------- volumes --------------------------

  // --- experimental hall (world volume)    ***LOOKHERE***
  //     beam line along the Z-axis
  G4double expHall_x = 10.0*m;  // half dimension along x 
  G4double expHall_y = 10.0*m;  // half dimension along y
  G4double expHall_z = 10.0*m;  // half dimension along z

  G4Box* experimentalHall_box
    = new G4Box( "expHall_box", expHall_x, expHall_y, expHall_z );

  experimentalHall_log = new G4LogicalVolume( experimentalHall_box, // solid 
                                              Vacuum,               // material
                                              "expHall_log",        // name
                                              0,                    // field manager
                                              0,                    // sensitive detector
                                              0 );                  // user limits

  experimentalHall_phys = new G4PVPlacement( 0,                     // rotation
                                             G4ThreeVector(),       // translation
                                             "expHall",             // name
                                             experimentalHall_log,  // logical volume
                                             0,                     // mother phys volume
                                             false,                 // boolean operation
                                             0 );                   // copy number

  //***LOOKHERE***  
  // Select before at which  z  the detector starts (by default,
  // at z = 0.0 ) , and then it extends to larger z .
  // You can also select event (air) gaps between:
  //   Tracker - EM Cal ;  EM Cal - HAD Cal ;  HAD Cal - muon detector 
  const G4double zStartDetector = 0.0;
  const G4double zGapTrackerEmCal = 0.0;
  const G4double zGapEmCalHadCal = 0.0;
  const G4double zGapHadCalMuonDetector = 0.0;

  G4double zPos = zStartDetector;

  G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
         << G4endl << "\t Before the Tracker " << G4endl
         << "\t zPos = " << zPos / m << " m " << G4endl; //***DEBUG***

  // --- First subdetector: the Tracker ---
  if ( theTrackerLength > 1.0E-06*mm ) {
    G4Tubs* solidTracker = new G4Tubs( "solidTracker", 
				       0.0,                  // inner radius
				       theDetectorRadius,    // outer radius
				       theTrackerLength/2.0, // half cylinder length in z
				       0.0,                  // starting phi angle in rad
				       2.0*pi );             // final phi angle in rad
    
    logicTracker = new G4LogicalVolume( solidTracker,        // solid 
					theTrackerMaterial,  // material
					"logicTracker",      // name
					0,                   // field manager
					0,                   // sensitive detector
					0 );                 // user limits
    
    zPos += theTrackerLength/2.0;
    physiTracker = new G4PVPlacement( 0,                       // rotation
				      G4ThreeVector(0.0,0.0,zPos), // translation
				      logicTracker,            // logical volume
				      "physiTracker",          // name
				      experimentalHall_log,    // mother logical volume
				      false,                   // boolean operation
				      0 );                     // copy number

    zPos += theTrackerLength/2.0;
  }

  G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
	 << G4endl << "\t After the Tracker " << G4endl
	 << "\t zPos = " << zPos / m << " m " << G4endl; //***DEBUG***

  zPos += zGapTrackerEmCal;

  // --- Second subdetector: the EM Calorimeter ---
  if ( theEmAbsorberTotalLength > 1.0E-06*mm  &&  theIsEmCalHomogeneous ) { 
    G4Tubs* solidEmCalo = new G4Tubs( "solidEmCalo", 
				       0.0,                  // inner radius
				       theDetectorRadius,    // outer radius
				       theEmAbsorberTotalLength/2.0, // half cylinder
				       0.0,                  // starting phi angle in rad
				       2.0*pi );             // final phi angle in rad
    
    logicEmCalo = new G4LogicalVolume( solidEmCalo  ,          // solid 
					theEmAbsorberMaterial, // material
					"logicEmCalo",         // name
					0,                     // field manager
					0,                     // sensitive detector
					0 );                   // user limits
    
    zPos += theEmAbsorberTotalLength/2.0;
    physiEmCalo = new G4PVPlacement( 0,                       // rotation
				     G4ThreeVector(0.0,0.0,zPos), // translation
				     logicEmCalo,             // logical volume
				     "physiEmCalo",           // name
				     experimentalHall_log,    // mother logical volume
				     false,                   // boolean operation
				     0 );                     // copy number

    zPos += theEmAbsorberTotalLength/2.0;
    
  } else if ( theEmAbsorberTotalLength > 1.0E-06*mm ) { 

    // For a sampling calorimeter, we use Replica placement. 
    // To do that, we have to define two extra volumes: the "calorimeter" volume 
    // and the "module". The former, which has the world as its mother volume, 
    // is the mother of the module volume. The calorimeter volume is completely 
    // filled by a number (theEmActiveLayerNumber) of replicas of the module volume. 
    // A module volume, in its turn, is the mother volume of the absorber layer + 
    // active layer. 

    // Absorber layer : logical
    G4double zAbsorber = theEmAbsorberTotalLength / 
      static_cast< double >( theEmActiveLayerNumber );
    zAbsorber /= 2.0;  // Half dimension along z .
    
    G4Tubs* solidEmAbsorber = new G4Tubs( "solidEmAbsorber", 
					  0.0,                // inner radius
					  theDetectorRadius,  // outer radius
					  zAbsorber,          // half cylinder length
					  0.0,                // starting phi angle 
					  2.0*pi );           // final phi angle in rad
    
    logicEmAbsorber = new G4LogicalVolume( solidEmAbsorber,       // solid 
					   theEmAbsorberMaterial, // material
					   "logicEmAbsorber",     // name
					   0,                     // field manager
					   0,                     // sensitive detector
					   0 );                   // user limits
    
    // Active layer : logical
    G4double zActive = theEmActiveLayerSize / 2.0;  // half dimension along z 
    
    G4Tubs* solidEmActive = new G4Tubs( "solidEmActive", 
					0.0,                // inner radius
					theDetectorRadius,  // outer radius
					zActive,            // half cylinder length in z
					0.0,                // starting phi angle in rad
					2.0*pi );           // final phi angle in rad
    
    logicEmActive = new G4LogicalVolume( solidEmActive,         // solid 
					 theEmActiveMaterial,   // material
					 "logicEmActive",       // name
					 0,                     // field manager
					 0,                     // sensitive detector
					 0 );                   // user limits
    
    // Module : logical
    G4double zModule = zAbsorber + zActive;  // half dimension along z 
    
    G4Tubs* solidEmModule = new G4Tubs( "solidEmModule", 
					0.0,                // inner radius
					theDetectorRadius,  // outer radius
					zModule,            // half cylinder length in z
					0.0,                // starting phi angle in rad
					2.0*pi );           // final phi angle in rad
    
    logicEmModule = new G4LogicalVolume( solidEmModule,   // solid 
					 Lead,            // material, it does NOT matter
					 "logicEmModule", // name
					 0,               // field manager
					 0,               // sensitive detector
					 0 );             // user limits
    
    // EM alorimeter : logical
    G4double zCalo = theEmActiveLayerNumber*zModule;  // half dimension along z 
    
    G4Tubs* solidEmCalo = new G4Tubs( "solidEmCalo", 
				      0.0,                // inner radius
				      theDetectorRadius,  // outer radius
				      zCalo,              // half cylinder length in z
				      0.0,                // starting phi angle in rad
				      2.0*pi );           // final phi angle in rad
    
    logicEmCalo = new G4LogicalVolume( solidEmCalo,      // solid 
				       Lead,             // material, it does NOT matter
				       "logicEmCalo",    // name
				       0,                // field manager
				       0,                // sensitive detector
				       0 );              // user limits
    
    // Absorber layer : physical
    G4double zpos = - zActive;
    physiEmAbsorber = new G4PVPlacement( 0,                       // rotation
					 G4ThreeVector(0.0,0.0,zpos), // translation
					 logicEmAbsorber,         // logical volume
					 "physiEmAbsorber",       // name
					 logicEmModule,           // mother logic volume
					 false,                   // boolean operation
					 0 );                     // copy number
    
    // Active layer : physical
    zpos += zAbsorber + zActive;
    physiEmActive = new G4PVPlacement( 0,                       // rotation
				       G4ThreeVector(0.0,0.0,zpos), // translation
				       logicEmActive,           // logical volume
				       "physiEmActive",         // name
				       logicEmModule,           // mother logical volume
				       false,                   // boolean operation
				       0 );                     // copy number
    
    // Module : physical (using replica)
    physiEmModule = new G4PVReplica( "EMCalo",                // name
				     logicEmModule,           // logical volume
				     logicEmCalo,             // mother logical volume
				     kZAxis,                  // axis of replication
				     theEmActiveLayerNumber,  // number of replica
				     2*(zAbsorber+zActive) ); // (full) width of replica
    
    // Calorimeter : physical
    zPos += zCalo;
    physiEmCalo = new G4PVPlacement( 0,                       // rotation
				     G4ThreeVector(0.0,0.0,zPos), // translation
				     "physiEmCalo",           // its name
				     logicEmCalo,             // logical volume
				     experimentalHall_phys,   // mother physical volume
				     false,                   // boolean operation
				     0 );                     // copy number
    
    zPos += zCalo;

    G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
           << G4endl << "\t EM Calorimeter : " << G4endl
           << "\t zAbsorber       = " << zAbsorber / mm <<  " mm " << G4endl
	   << "\t zActive         = " << zActive / mm <<  " mm " << G4endl
	   << "\t zModule         = " << zModule / mm << " mm " << G4endl
	   << "\t total Absorber  = " << zAbsorber*theEmActiveLayerNumber / m 
           << " m " << G4endl
	   << "\t zCalo           = " << zCalo / m << " m " << G4endl; //***DEBUG***

  }

  G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
	 << G4endl << "\t After the EM Calorimeter " << G4endl
	 << "\t zPos = " << zPos / m << " m " << G4endl; //***DEBUG***

  zPos += zGapEmCalHadCal;

  // --- Third subdetector: the HAD Calorimeter ---
  if ( theHadAbsorberTotalLength > 1.0E-06*mm  &&  theIsHadCalHomogeneous ) { 
    G4Tubs* solidHadCalo = new G4Tubs( "solidHadCalo", 
				       0.0,                  // inner radius
				       theDetectorRadius,    // outer radius
				       theHadAbsorberTotalLength/2.0, // half cylinder
				       0.0,                  // starting phi angle in rad
				       2.0*pi );             // final phi angle in rad
    
    logicHadCalo = new G4LogicalVolume( solidHadCalo  ,         // solid 
					theHadAbsorberMaterial, // material
					"logicHadCalo",         // name
					0,                      // field manager
					0,                      // sensitive detector
					0 );                    // user limits
    
    zPos += theHadAbsorberTotalLength/2.0;
    physiHadCalo = new G4PVPlacement( 0,                         // rotation
				      G4ThreeVector(0.0,0.0,zPos), // translation
				      logicHadCalo,              // logical volume
				      "physiHadCalo",            // name
				      experimentalHall_log,      // mother logical volume
				      false,                     // boolean operation
				      0 );                       // copy number

    zPos += theHadAbsorberTotalLength/2.0;
    
  } else if ( theHadAbsorberTotalLength > 1.0E-06*mm ) { 
 
    // We use Replica as for the Electromagnetic calorimeter.

    // Absorber layer : logical
    G4double zAbsorber = theHadAbsorberTotalLength / 
      static_cast< double >( theHadActiveLayerNumber );
    zAbsorber /= 2.0;  // Half dimension along z .
    
    G4Tubs* solidHadAbsorber = new G4Tubs( "solidHadAbsorber", 
					   0.0,                // inner radius
					   theDetectorRadius,  // outer radius
					   zAbsorber,          // half cylinder length
					   0.0,                // starting phi angle 
					   2.0*pi );           // final phi angle in rad
    
    logicHadAbsorber = new G4LogicalVolume( solidHadAbsorber,       // solid 
					    theHadAbsorberMaterial, // material
					    "logicHadAbsorber",     // name
					    0,                      // field manager
					    0,                      // sensitive detector
					    0 );                    // user limits
    
    // Active layer : logical
    G4double zActive = theHadActiveLayerSize / 2.0;  // half dimension along z 
    
    G4Tubs* solidHadActive = new G4Tubs( "solidHadActive", 
					 0.0,                // inner radius
					 theDetectorRadius,  // outer radius
					 zActive,            // half cylinder length in z
					 0.0,                // starting phi angle in rad
					 2.0*pi );           // final phi angle in rad
    
    logicHadActive = new G4LogicalVolume( solidHadActive,         // solid 
					  theHadActiveMaterial,   // material
					  "logicHadActive",       // name
					  0,                      // field manager
					  0,                      // sensitive detector
					  0 );                    // user limits
    
    // Module : logical
    G4double zModule = zAbsorber + zActive;  // half dimension along z 
    
    G4Tubs* solidHadModule = new G4Tubs( "solidHadModule", 
					 0.0,                // inner radius
					 theDetectorRadius,  // outer radius
					 zModule,            // half cylinder length in z
					 0.0,                // starting phi angle in rad
					 2.0*pi );           // final phi angle in rad
    
    logicHadModule = new G4LogicalVolume( solidHadModule,   // solid 
					  Lead,           // material, it does NOT matter
					  "logicHadModule", // name
					  0,                // field manager
					  0,                // sensitive detector
					  0 );              // user limits
    
    // HAD calorimeter : logical
    G4double zCalo = theHadActiveLayerNumber*zModule;  // half dimension along z 
    
    G4Tubs* solidHadCalo = new G4Tubs( "solidHadCalo", 
				       0.0,                // inner radius
				       theDetectorRadius,  // outer radius
				       zCalo,              // half cylinder length in z
				       0.0,                // starting phi angle in rad
				       2.0*pi );           // final phi angle in rad
    
    logicHadCalo = new G4LogicalVolume( solidHadCalo,     // solid 
				        Lead,             // material, it does NOT matter
				        "logicHadCalo",   // name
				        0,                // field manager
				        0,                // sensitive detector
				        0 );              // user limits
    
    // Absorber layer : physical
    G4double zpos = - zActive;
    physiHadAbsorber = new G4PVPlacement( 0,                       // rotation
					  G4ThreeVector(0.0,0.0,zpos), // translation
					  logicHadAbsorber,        // logical volume
					  "physiHadAbsorber",      // name
					  logicHadModule,          // mother logic volume
					  false,                   // boolean operation
					  1000 );                  // copy number
    
    // Active layer : physical
    zpos += zAbsorber + zActive;
    physiHadActive = new G4PVPlacement( 0,                       // rotation
				        G4ThreeVector(0.0,0.0,zpos), // translation
				        logicHadActive,          // logical volume
				        "physiHadActive",        // name
				        logicHadModule,          // mother logical volume
				        false,                   // boolean operation
				        0 );                     // copy number
    
    // Module : physical (using replica)
    physiHadModule = new G4PVReplica( "HADCalo",               // name
				      logicHadModule,          // logical volume
				      logicHadCalo,            // mother logical volume
				      kZAxis,                  // axis of replication
				      theHadActiveLayerNumber, // number of replica
				      2*(zAbsorber+zActive) ); // (full) width of replica
    
    // Calorimeter : physical
    zPos += zCalo;
    physiHadCalo = new G4PVPlacement( 0,                       // rotation
				      G4ThreeVector(0.0,0.0,zPos), // translation
				      "physiHadCalo",          // its name
				      logicHadCalo,            // logical volume
				      experimentalHall_phys,   // mother physical volume
				      false,                   // boolean operation
				      0 );                     // copy number
    
    zPos += zCalo;

    G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
           << G4endl << "\t HAD Calorimeter : " << G4endl
           << "\t zAbsorber      = " << zAbsorber / mm <<  " mm " << G4endl
           << "\t zActive        = " << zActive / mm <<  " mm " << G4endl
	   << "\t zModule        = " << zModule / mm << " mm " << G4endl
	   << "\t total Absorber = " << zAbsorber*theHadActiveLayerNumber / m 
           << " m " << G4endl
	   << "\t zCalo          = " << zCalo / m << " m " << G4endl; //***DEBUG***

  }

  G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
	 << G4endl << "\t After the HAD Calorimeter " << G4endl
	 << "\t zPos = " << zPos / m << " m " << G4endl; //***DEBUG***

  zPos += zGapHadCalMuonDetector;

  // --- Fourth (last) subdetector: the Muon detector ---
  if ( theMuonLength > 1.0E-06*mm ) {
    G4Tubs* solidMuon = new G4Tubs( "solidMuon", 
				    0.0,                 // inner radius
				    theDetectorRadius,   // outer radius
				    theMuonLength/2.0,   // half cylinder length in z
				    0.0,                 // starting phi angle in rad
				    2.0*pi );            // final phi angle in rad
    
    logicMuon = new G4LogicalVolume( solidMuon,          // solid 
				     theMuonMaterial,    // material
				     "logicMuon",        // name
				     0,                  // field manager
				     0,                  // sensitive detector
				     0 );                // user limits
    
    zPos += theMuonLength/2.0;
    physiMuon = new G4PVPlacement( 0,                       // rotation
				   G4ThreeVector(0.0,0.0,zPos), // translation
				   logicMuon,               // logical volume
				   "physiMuon",             // name
				   experimentalHall_log,    // mother logical volume
				   false,                   // boolean operation
				   0 );                     // copy number

    zPos += theMuonLength/2.0;
  }

  G4cout << " DetectorSliceDetectorConstruction::ConstructCalorimeter() : DEBUG Info "
	 << G4endl << "\t After the Muon Detector " << G4endl
	 << "\t zPos = " << zPos / m << " m " << G4endl; //***DEBUG***

  // --- Sensitive detectors

  if ( ! theSensitiveEmCalorimeter ) { 
    theSensitiveEmCalorimeter = 
      new DetectorSliceSensitiveEmCalo( "sensitiveEmCalorimeter" );
  }
  G4SDManager::GetSDMpointer()->AddNewDetector( theSensitiveEmCalorimeter );
  if ( theIsEmCalHomogeneous ) {
    // Trick to have the information on the visible energy per
    // particle type when the electromagnetic calorimeter is
    // homogeneous (e.g. PbWO4).
    if (  logicEmCalo ) {
      logicEmCalo->SetSensitiveDetector( theSensitiveEmCalorimeter );
    }
  } else if ( logicEmActive ) {
    logicEmActive->SetSensitiveDetector( theSensitiveEmCalorimeter );
  }

  if ( ! theSensitiveHadCalorimeter ) { 
    theSensitiveHadCalorimeter = 
      new DetectorSliceSensitiveHadCalo( "sensitiveHadCalorimeter" );
  }
  G4SDManager::GetSDMpointer()->AddNewDetector( theSensitiveHadCalorimeter );
  if ( theIsHadCalHomogeneous ) {
    if ( logicHadCalo ) {
      logicHadCalo->SetSensitiveDetector( theSensitiveHadCalorimeter );
    }
  } else if ( logicHadActive ) {
    logicHadActive->SetSensitiveDetector( theSensitiveHadCalorimeter );
  }

  // --- Visualization attributes
  experimentalHall_log->SetVisAttributes( G4VisAttributes::Invisible );
  // The World is not visualized.

  logicEmCalo->SetVisAttributes( G4VisAttributes::Invisible );
  logicHadCalo->SetVisAttributes( G4VisAttributes::Invisible );
  // The calo is not visualized.

  logicEmModule->SetVisAttributes( G4VisAttributes::Invisible );
  logicHadModule->SetVisAttributes( G4VisAttributes::Invisible );
  // The module is not visualized.

  if ( ! theVisAttAbsorber ) {
    theVisAttAbsorber = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
    theVisAttAbsorber->SetVisibility( true );
    theVisAttAbsorber->SetForceWireframe( true );
  }
  logicEmAbsorber->SetVisAttributes( theVisAttAbsorber );
  logicHadAbsorber->SetVisAttributes( theVisAttAbsorber );
  // The absorber layer will appear in white colour.
  // (the order of colours is: (red, green, blue) )
  
  if ( ! theVisAttActive ) {
    theVisAttActive = new G4VisAttributes( G4Colour( 1.0, 1.0, 0.0 ) );
    theVisAttActive->SetVisibility( true );
    theVisAttActive->SetForceWireframe( true );
  }
  logicEmActive->SetVisAttributes( theVisAttActive );
  logicHadActive->SetVisAttributes( theVisAttActive );
  // The active layer will appear in yellow colour.
  
  //G4cout << " END  DetectorSliceDetectorConstruction::ConstructCalorimeter()" 
  //       << G4endl; //***DEBUG***

  return experimentalHall_phys;
}


G4bool DetectorSliceDetectorConstruction::areParametersOK() {

  bool isOk = true;

  if ( ! theTrackerMaterial ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : UNDEFINED tracker material" << G4endl;
  }
  if ( theTrackerLength < 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theTrackerLength = " << theTrackerLength << G4endl;
  }

  if ( ! theEmAbsorberMaterial ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : UNDEFINED EM absorber material" << G4endl;
  }
  if ( ! theEmActiveMaterial ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : UNDEFINED EM active material" << G4endl;
  }
  if ( theEmAbsorberTotalLength < 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theEmAbsorberTotalLength = " << theEmAbsorberTotalLength << G4endl;
  }
  if ( theEmActiveLayerNumber < 0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theEmActiveLayerNumber = " << theEmActiveLayerNumber << G4endl;
  }
  if ( theEmActiveLayerSize < 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theEmActiveLayerSize = " << theEmActiveLayerSize << G4endl;
  }

  if ( ! theHadAbsorberMaterial ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : UNDEFINED HAD absorber material" << G4endl;
  }
  if ( ! theHadActiveMaterial ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : UNDEFINED HAD active material" << G4endl;
  }
  if ( theHadAbsorberTotalLength < 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theHadAbsorberTotalLength = " << theHadAbsorberTotalLength << G4endl;
  }
  if ( theHadActiveLayerNumber < 0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theHadActiveLayerNumber = " << theHadActiveLayerNumber << G4endl;
  }
  if ( theHadActiveLayerSize < 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theHadActiveLayerSize = " << theHadActiveLayerSize << G4endl;
  }

  if ( ! theMuonMaterial ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : UNDEFINED muon material" << G4endl;
  }
  if ( theMuonLength < 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theMuonLength = " << theMuonLength << G4endl;
  }

  if ( theDetectorRadius <= 0.0 ) {
    isOk = false;
    G4cout << " DetectorSliceDetectorConstruction::areParametersOK() : theDetectorRadius = " << theDetectorRadius << G4endl;
  }

  return isOk;
}


void DetectorSliceDetectorConstruction::SetTrackerMaterial( const G4String name ) {

  if ( name == "Scintillator" || name == "scintillator" ) { 
    theTrackerMaterial = Polystyrene;
  } else if ( name == "LAr" ||
 	      name == "LiquidArgon" || name == "liquidArgon" ) { 
    theTrackerMaterial = LiquidArgon;
  } else if ( name == "Si" ||
 	      name == "Silicon" || name == "silicon" ) { 
    theTrackerMaterial = Silicon;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Silicon *  will be used." 
	   << G4endl << G4endl;  
    theTrackerMaterial = Silicon;
  }
  
  logicTracker->SetMaterial( theTrackerMaterial );
  
  G4cout << " Tracker Material = " << logicTracker->GetMaterial()->GetName() 
         << G4endl;
    
}


void DetectorSliceDetectorConstruction::SetEmAbsorberMaterial( const G4String name ) {

  if ( name == "Fe" ||
       name == "Iron" || name == "iron" ) { 
    theEmAbsorberMaterial = Iron;
  } else if ( name == "Cu" ||
 	      name == "Copper" || name == "copper" ) { 
    theEmAbsorberMaterial = Copper;
  } else if ( name == "Pb" ||
 	      name == "Lead" || name == "lead" ) { 
    theEmAbsorberMaterial = Lead;
  } else if ( name == "PbWO4" ) {
    theEmAbsorberMaterial = PbWO4;
  } else if ( name == "W" ||
 	      name == "Tungsten" || name == "tungsten" ) { 
    theEmAbsorberMaterial = Tungsten;
  } else if ( name == "U" ||
 	      name == "Uranium" || name == "uranium" ) { 
    theEmAbsorberMaterial = Uranium;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Lead *  will be used." 
	   << G4endl << G4endl;  
    theEmAbsorberMaterial = Lead;
  }
  
  logicEmAbsorber->SetMaterial( theEmAbsorberMaterial );
  
  G4cout << " EM Absorber Material = " << logicEmAbsorber->GetMaterial()->GetName() 
         << G4endl;
  
}


void DetectorSliceDetectorConstruction::SetEmActiveMaterial( const G4String name ) {

  if ( name == "Scintillator" || name == "scintillator" ) { 
    theEmActiveMaterial = Polystyrene;
  } else if ( name == "LAr" ||
 	      name == "LiquidArgon" || name == "liquidArgon" ) { 
    theEmActiveMaterial = LiquidArgon;
  } else if ( name == "PbWO4" ) {
    theEmActiveMaterial = PbWO4;
  } else if ( name == "Si" ||
 	      name == "Silicon" || name == "silicon" ) { 
    theEmActiveMaterial = Silicon;
  } else if ( name == "Quartz" || name == "quartz" ) { 
    theEmActiveMaterial = Quartz;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * LiquidArgon *  will be used." 
	   << G4endl << G4endl;  
    theEmActiveMaterial = LiquidArgon;
  }
  
  logicEmActive->SetMaterial( theEmActiveMaterial );
  
  G4cout << " EM Active Material = " << logicEmActive->GetMaterial()->GetName() 
         << G4endl;
  
}


void DetectorSliceDetectorConstruction::SetHadAbsorberMaterial( const G4String name ) {

  if ( name == "Fe" ||
       name == "Iron" || name == "iron" ) { 
    theHadAbsorberMaterial = Iron;
  } else if ( name == "Cu" ||
 	      name == "Copper" || name == "copper" ) { 
    theHadAbsorberMaterial = Copper;
  } else if ( name == "Pb" ||
 	      name == "Lead" || name == "lead" ) { 
    theHadAbsorberMaterial = Lead;
  } else if ( name == "PbWO4" ) {
    theHadAbsorberMaterial = PbWO4;
  } else if ( name == "W" ||
 	      name == "Tungsten" || name == "tungsten" ) { 
    theHadAbsorberMaterial = Tungsten;
  } else if ( name == "U" ||
 	      name == "Uranium" || name == "uranium" ) { 
    theHadAbsorberMaterial = Uranium;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Iron *  will be used." 
	   << G4endl << G4endl;  
    theHadAbsorberMaterial = Iron;
  }
  
  logicHadAbsorber->SetMaterial( theHadAbsorberMaterial );
  
  G4cout << " HAD Absorber Material = " << logicHadAbsorber->GetMaterial()->GetName() 
         << G4endl;
  
}


void DetectorSliceDetectorConstruction::SetHadActiveMaterial( const G4String name ) {

  if ( name == "Scintillator" || name == "scintillator" ) { 
    theHadActiveMaterial = Polystyrene;
  } else if ( name == "LAr" ||
 	      name == "LiquidArgon" || name == "liquidArgon" ) { 
    theHadActiveMaterial = LiquidArgon;
  } else if ( name == "PbWO4" ) {
    theHadActiveMaterial = PbWO4;
  } else if ( name == "Si" ||
 	      name == "Silicon" || name == "silicon" ) { 
    theHadActiveMaterial = Silicon;
  } else if ( name == "Quartz" || name == "quartz" ) { 
    theHadActiveMaterial = Quartz;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Scintillator *  will be used." 
	   << G4endl << G4endl;  
    theHadActiveMaterial = Polystyrene;
  }
  
  logicHadActive->SetMaterial( theHadActiveMaterial );
  
  G4cout << " HAD Active Material = " << logicHadActive->GetMaterial()->GetName() 
         << G4endl;
  
}


void DetectorSliceDetectorConstruction::SetMuonMaterial( const G4String name ) {

  if ( name == "Fe" ||
       name == "Iron" || name == "iron" ) { 
    theMuonMaterial = Iron;
  } else if ( name == "Cu" ||
 	      name == "Copper" || name == "copper" ) { 
    theMuonMaterial = Copper;
  } else if ( name == "Pb" ||
 	      name == "Lead" || name == "lead" ) { 
    theMuonMaterial = Lead;
  } else if ( name == "W" ||
 	      name == "Tungsten" || name == "tungsten" ) { 
    theMuonMaterial = Tungsten;
  } else if ( name == "U" ||
 	      name == "Uranium" || name == "uranium" ) { 
    theHadAbsorberMaterial = Uranium;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Iron *  will be used." 
	   << G4endl << G4endl;  
    theHadAbsorberMaterial = Iron;
  }
    
  logicMuon->SetMaterial( theMuonMaterial );
  
  G4cout << " Muon Material = " << logicMuon->GetMaterial()->GetName() 
         << G4endl;
    
}


void DetectorSliceDetectorConstruction::UpdateGeometry() {

  //G4cout << " BEGIN  DetectorSliceDetectorConstruction::UpdateGeometry" << G4endl; //***DEBUG***

  G4RunManager::GetRunManager()->DefineWorldVolume( ConstructCalorimeter() );

  PrintParameters();

  //G4cout << " END  DetectorSliceDetectorConstruction::UpdateGeometry" << G4endl; //***DEBUG***
}


void DetectorSliceDetectorConstruction::PrintParameters() {

  G4cout << G4endl << G4endl
         << " ------  DetectorSliceDetectorConstruction::PrintParameters() ------ " 
	 << G4endl;

  G4cout << " Tracker : " << G4endl 
         << "\t material : ";
  if ( theTrackerMaterial ) {
    G4cout << theTrackerMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << "\t Tracker Length = "
	 << theTrackerLength / mm << " mm" << G4endl;

  G4cout << " EM Calorimeter : " << G4endl 
         << "\t absorber : ";
  if ( theEmAbsorberMaterial ) {
    G4cout << theEmAbsorberMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << "\t active Material: ";
  if ( theEmActiveMaterial ) {
    G4cout << theEmActiveMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << "\t Is Homogeneous ? " << theIsEmCalHomogeneous;
  G4cout << G4endl << "\t Absorber Total Length = "
	 << theEmAbsorberTotalLength / mm << " mm";
  G4cout << G4endl << "\t Active Layer Number   = " << theEmActiveLayerNumber;
  G4cout << G4endl << "\t Active Layer Size     = " << theEmActiveLayerSize/mm 
	 << " mm" << G4endl;

  G4cout << " HAD Calorimeter : " << G4endl 
         << "\t absorber : ";
  if ( theHadAbsorberMaterial ) {
    G4cout << theHadAbsorberMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << "\t active Material: ";
  if ( theHadActiveMaterial ) {
    G4cout << theHadActiveMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << "\t Is Homogeneous ? " << theIsHadCalHomogeneous;
  G4cout << G4endl << "\t Absorber Total Length = "
	 << theHadAbsorberTotalLength / mm << " mm";
  G4cout << G4endl << "\t Active Layer Number   = " << theHadActiveLayerNumber;
  G4cout << G4endl << "\t Active Layer Size     = " << theHadActiveLayerSize/mm 
	 << " mm" << G4endl;

  G4cout << " Muon detector : " << G4endl 
         << "\t material : ";
  if ( theMuonMaterial ) {
    G4cout << theMuonMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << "\t Muon detector Length = "
	 << theMuonLength / mm << " mm" << G4endl;

  G4cout << " Detector Radius = " << theDetectorRadius / mm << " mm";

  G4cout << G4endl << " -------------------------------------------------------- "
         << G4endl << G4endl;

}
