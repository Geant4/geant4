#include "MyDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
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

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
      
#include "MyDetectorMessenger.hh"

#include "MySensitiveCalorimeter.hh"
#include "G4SDManager.hh"
            

MyDetectorConstruction::MyDetectorConstruction() :  
  Iron(0), Copper(0), Tungsten(0), Lead(0), Uranium(0), PbWO4(0), 
  Polystyrene(0), LiquidArgon(0), Silicon(0), Quartz(0),
  theAbsorberMaterial(0), theActiveMaterial(0),
  experimentalHall_log(0), experimentalHall_phys(0),
  logicCalo(0), physiCalo(0),
  logicModule(0), physiModule(0),
  logicAbsorber(0), physiAbsorber(0),
  logicActive(0), physiActive(0),
  uniformMagField(0)
{
  G4FieldManager* fieldMgr 
    = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  if ( uniformMagField ) delete uniformMagField;

  // Apply a global uniform magnetic field along Z axis.
  uniformMagField = new G4UniformMagField( G4ThreeVector(0.0, 0.0, 0.0*tesla) );
  fieldMgr->SetDetectorField( uniformMagField );
  fieldMgr->CreateChordFinder( uniformMagField );
  
  detectorMessenger = new MyDetectorMessenger(this);
}


MyDetectorConstruction::~MyDetectorConstruction() {
  delete detectorMessenger;
  delete uniformMagField;
}


G4VPhysicalVolume* MyDetectorConstruction::Construct() {

  //------------------- materials ------------------------

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, pressure, temperature, fractionmass;
  G4String name, symbol;
  G4int nel, natoms;

  //--- elements

  a = 1.01*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H2", z=1., a);

  a = 2.01*g/mole;
  G4Element* elD = new G4Element(name="Deuterium", symbol="D", z=1., a);

  a = 4.*g/mole;
  G4Element* elHe = new G4Element(name="Helium", symbol="He", z=2., a);

  a = 6.94*g/mole;
  G4Element* elLi = new G4Element(name="Lithium", symbol="Li", z=3., a);

  a = 9.01*g/mole;
  G4Element* elBe = new G4Element(name="Berillium", symbol="Be", z=4., a);

  a = 12.01*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a);

  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N2", z=7., a);

  a = 16.*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O2", z=8., a);

  a = 20.18*g/mole;
  G4Element* elNe = new G4Element(name="Neon", symbol="Ne", z=10., a);

  a = 22.99*g/mole;
  G4Element* elNa = new G4Element(name="Sodium", symbol="Na", z=11., a);

  a = 26.98*g/mole;
  G4Element* elAl = new G4Element(name="Aluminium", symbol="Al", z=13., a);

  a = 28.085*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", z=14., a);

  a = 40.08*g/mole;
  G4Element* elCa = new G4Element(name="Calcium", symbol="Ca", z=20., a);

  a = 55.850*g/mole;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", z=26., a);

  a = 63.54*g/mole;
  G4Element* elCu = new G4Element(name="Copper", symbol="Cu", z=29., a);

  a = 183.85*g/mole;
  G4Element* elW = new G4Element(name="Tungstenm", symbol="W", z=74., a);

  a = 207.19*g/mole;
  G4Element* elPb = new G4Element(name="Lead", symbol="Pb", z=82., a);

  a = 238.03*g/mole;
  G4Element* elU = new G4Element(name="Uranium", symbol="U", z=92., a);

  //--- simple materials

  density = 2.7*g/cm3;
  a = 26.98*g/mole;
  G4Material* Aluminium = new G4Material(name="Aluminium", z=13., a, density);
  
  // Iron has a  X0 = 1.7585 cm  and  lambda_I = 16.760 cm.   
  density = 7.87*g/cm3;
  a = 55.85*g/mole;
  Iron = new G4Material(name="Iron", z=26., a, density);

  // Copper has a  X0 = 1.4353 cm  and  lambda_I = 15.056 cm.   
  density = 8.96*g/cm3;
  a = 63.54*g/mole;
  Copper = new G4Material(name="Copper", z=29., a, density);

  // Tungsten has a  X0 = 0.35 cm  and  lambda_I = 9.5855 cm. 
  density = 19.3*g/cm3;
  a = 183.85*g/mole;
  Tungsten = new G4Material(name="Tungsten", z=74., a, density);

  // Lead has a  X0 = 0.56120 cm  and  lambda_I = 17.092 cm.  
  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  Lead = new G4Material(name="Lead", z=82., a, density);

  // Uranium has a  X0 = 0.31662 cm  and  lambda_I = 10.501 cm.  
  density =  18.95*g/cm3;
  a = 238.03*g/mole;
  Uranium = new G4Material(name="Uranium", z=92., a, density);

  density = 1.4*g/cm3;
  a = 39.95*g/mole;
  LiquidArgon = new G4Material(name="LiquidArgon", z=18., a, density);

  density = 0.002*g/cm3;
  a = 39.95*g/mole;
  G4Material* ArgonGas = new G4Material(name="ArgonGas", z=18., a, density);

  density = 2.33*g/cm3;
  a = 28.085*g/mole;
  Silicon = new G4Material(name="Silicon", z=14., a, density);
  
  density = 8.96*g/cm3;
  a = 58.69*g/mole;
  G4Material* Nickel = new G4Material(name="Nickel", z=28., a, density);

  //--- mixtures

  density = 1.290*mg/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, 0.7);
  Air->AddElement(elO, 0.3);

  density     = 1.e-5*g/cm3;
  pressure    = 2.e-2*bar;
  temperature = STP_Temperature;  // From PhysicalConstants.h .
  G4Material* Vacuum = new G4Material(name="Vacuum", density, nel=1,
				      kStateGas, temperature, pressure);
  Vacuum->AddMaterial(Air, fractionmass=1.);

  // Plastic scintillator tiles (used both in CMS hadron calorimeter
  // and ATLAS hadron barrel calorimeter).
  density = 1.032*g/cm3;
  Polystyrene = new G4Material(name="Polystyrene", density, nel=2);
  Polystyrene->AddElement(elC, natoms=19);
  Polystyrene->AddElement(elH, natoms=21);

  // PbWO4 CMS crystals. It has a  X0 = 0.89 cm  and  lambda_I = 22.4 cm. 
  density = 8.28*g/cm3;
  PbWO4 = new G4Material(name="PbWO4", density, nel=3);
  PbWO4->AddElement(elPb, natoms=1);
  PbWO4->AddElement(elW,  natoms=1);
  PbWO4->AddElement(elO,  natoms=4);

  Quartz = new G4Material(name="Quartz", density=2.200*g/cm3, nel=2);
  Quartz->AddElement(elSi, 1);
  Quartz->AddElement(elO , 2);

  // --- Set default values
  theAbsorberMaterial = Iron;
  theActiveMaterial = Polystyrene;

  //------------------- volumes --------------------------

  // --- experimental hall (world volume)
  //     beam line along x axis
  G4double expHall_x = 2.0*m;  // half dimension along x 
  G4double expHall_y = 2.0*m;  // half dimension along y
  G4double expHall_z = 2.0*m;  // half dimension along z

  G4Box* experimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  experimentalHall_log = new G4LogicalVolume(experimentalHall_box, // solid 
                                             Vacuum,               // material
                                             "expHall_log",        // name
                                             0,                    // field manager
                                             0,                    // sensitive detector
                                             0);                   // user limits

  experimentalHall_phys = new G4PVPlacement(0,                     // rotation
                                            G4ThreeVector(),       // translation
                                            "expHall",             // name
                                            experimentalHall_log,  // logical volume
                                            0,                     // mother physical volume
                                            false,                 // boolean operation
                                            0);                    // copy number
  
  // --- Detector
  // The idea is to use Replica placement. 
  // To do that, we have to define two extra volumes: the "calorimeter" volume 
  // and the "module". The former, which has the world as its mother volume, 
  // is the mother of the module volume. The calorimeter volume is completely 
  // filled by 50 replicas of the module volume. A module volume, in its turn, 
  // is the mother volume of the absorber layer + active layer. 
  // We have chosen the following dimension of the sampling calorimeter:
  //   transversal dimensions:  2m x 2m   very large to have full lateral containment;  
  //   absorber thickness: 4 cm   similar to ATLAS and CMS tile hadronic calorimeters;
  //   active layer:       4 mm     "      "   "    "   "    "     "         "
  //   number of modules:  50  very large to have full longitudinal containement:
  //                           it corresponds to more than 10 interaction lengths
  //                           for: tungsten (20.9), uranium (19), 
  //                                copper (13), iron (12.5), lead (11.8), 
  //                                and about 9 for PbWO4 (which indeed
  //                                becomes 9.9 if we use PbWO4 also as active, as
  //                                fake for a real homogeneous calorimeter made 
  //                                entirely of PbWO4).

  //            --- absorber layer : logical
  G4double xAbsorber =   2.0*cm;  // half dimension along x 
  G4double yAbsorber = 100.0*cm;  // half dimension along y
  G4double zAbsorber = 100.0*cm;  // half dimension along z

  G4Box* solidAbsorber = new G4Box("solidAbsorber", xAbsorber, yAbsorber, zAbsorber);

  logicAbsorber = new G4LogicalVolume(solidAbsorber,       // solid 
                                      theAbsorberMaterial, // material
                                      "logicAbsorber",     // name
                                      0,                   // field manager
                                      0,                   // sensitive detector
                                      0);                  // user limits

  //            --- active layer : logical
  G4double xActive = 0.2*cm;     // half dimension along x 
  G4double yActive = yAbsorber;  // half dimension along y 
  G4double zActive = zAbsorber;  // half dimension along z 

  G4Box* solidActive = new G4Box("solidActive", xActive, yActive, zActive);

  logicActive = new G4LogicalVolume(solidActive,           // solid 
                                    theActiveMaterial,     // material
                                    "logicActive",         // name
                                    0,                     // field manager
                                    0,                     // sensitive detector
                                    0);                    // user limits

  //        --- module : logical
  G4double xModule = (xAbsorber+xActive);  // half dimension along x 
  G4double yModule = yAbsorber;            // half dimension along y
  G4double zModule = zAbsorber;            // half dimension along z

  G4Box* solidModule = new G4Box("solidModule", xModule, yModule, zModule);

  logicModule = new G4LogicalVolume(solidModule,    // solid 
                                    Lead,           // material, it does NOT matter
                                    "logicModule",  // name
                                    0,              // field manager
                                    0,              // sensitive detector
                                    0);             // user limits

  //    --- calorimeter : logical
  G4int numberOfModules = 50;
  G4double xCalo = numberOfModules*xModule;  // half dimension along x 
  G4double yCalo = yModule;                  // half dimension along y
  G4double zCalo = zModule;                  // half dimension along z

  G4Box* solidCalo = new G4Box("solidCalo", xCalo, yCalo, zCalo);

  logicCalo = new G4LogicalVolume(solidCalo,        // solid 
                                  Lead,             // material, it does NOT matter
                                  "logicCalo",      // name
                                  0,                // field manager
                                  0,                // sensitive detector
                                  0);               // user limits

  //            --- absorber layer : physical
  G4double xpos = - xActive;
  physiAbsorber = new G4PVPlacement(0,                       // rotation
                                    G4ThreeVector(xpos,0,0), // translation
                                    logicAbsorber,           // logical volume
                                    "physiAbsorber",         // name
                                    logicModule,             // mother logical volume
                                    false,                   // boolean operation
                                    1000);                   // copy number

  //            --- active layer : physical
  xpos += xAbsorber + xActive;
  physiActive = new G4PVPlacement(0,                       // rotation
                                  G4ThreeVector(xpos,0,0), // translation
                                  logicActive,             // logical volume
                                  "physiActive",           // name
                                  logicModule,             // mother logical volume
                                  false,                   // boolean operation
                                  2000);                   // copy number

  //        --- module : physical (using replica)
  physiModule = new G4PVReplica("Calo",                  // name
                                logicModule,             // logical volume
                                logicCalo,               // mother logical volume
                                kXAxis,                  // axis of replication
                                numberOfModules,         // number of replica
                                2*(xAbsorber+xActive) ); // (full) width of replica

  //    --- calorimeter : physical
  physiCalo = new G4PVPlacement(0,                     // rotation
                                G4ThreeVector(),       // translation
                                "physiCalo",           // its name
                                logicCalo,             // logical volume
                                experimentalHall_phys, // mother physical volume
                                false,                 // boolean operation
                                100);                  // copy number

  // --- Sensitive detectors
  G4String sensitiveCalorimeterName = "My/SensitiveCalorimeter";
  MySensitiveCalorimeter* aSensitiveCalorimeter = 
    new MySensitiveCalorimeter( sensitiveCalorimeterName );
  G4SDManager::GetSDMpointer()->AddNewDetector( aSensitiveCalorimeter );
  logicActive->SetSensitiveDetector( aSensitiveCalorimeter );
  
  // --- Visualization attributes
  experimentalHall_log->SetVisAttributes( G4VisAttributes::Invisible );
  // The World is not visualized.

  logicCalo->SetVisAttributes( G4VisAttributes::Invisible );
  // The calo is not visualized.

  logicModule->SetVisAttributes( G4VisAttributes::Invisible );
  // The calo is not visualized.

  G4VisAttributes* visAttAbsorber = new G4VisAttributes( G4Colour(1.0,1.0,1.0) );
  visAttAbsorber->SetVisibility(true);
  visAttAbsorber->SetForceWireframe(true);
  logicAbsorber->SetVisAttributes(visAttAbsorber);
  // The absorber layer will appear in white colour.
  // (the order of colours is: (red, green, blue) )
  
  G4VisAttributes* visAttActive = new G4VisAttributes( G4Colour(1.0,1.0,0.0) );
  visAttActive->SetVisibility(true);
  visAttActive->SetForceWireframe(true);
  logicActive->SetVisAttributes(visAttActive);
  // The active layer will appear in yellow colour.
  
  //------------------------------------------------------

  PrintParameters();

  return experimentalHall_phys;
}


void MyDetectorConstruction::PrintParameters() {

  G4cout << G4endl << G4endl
         << " ------  MyDetectorConstruction::PrintParameters() ------ " << G4endl
         << " Absorber Material = ";
  if ( theAbsorberMaterial ) {
    G4cout << theAbsorberMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << " Active Material   = ";
  if ( theActiveMaterial ) {
    G4cout << theActiveMaterial->GetName();
  } else {
    G4cout << " UNDEFINED ";
  }
  G4cout << G4endl << " -------------------------------------------------------- "
         << G4endl << G4endl;

}


void MyDetectorConstruction::SetMagField(G4double fieldValue) {
  uniformMagField->SetFieldValue( G4ThreeVector(0.0, 0.0, fieldValue) );
}


void MyDetectorConstruction::SetAbsorberMaterial(const G4String name) {

  if ( name == "Fe" ||
       name == "Iron" || name == "iron" ) { 
    theAbsorberMaterial = Iron;
  } else if ( name == "Cu" ||
 	      name == "Copper" || name == "copper" ) { 
    theAbsorberMaterial = Copper;
  } else if ( name == "Pb" ||
 	      name == "Lead" || name == "lead" ) { 
    theAbsorberMaterial = Lead;
  } else if ( name == "PbWO4" ) {
    theAbsorberMaterial = PbWO4;
  } else if ( name == "W" ||
 	      name == "Tungsten" || name == "tungsten" ) { 
    theAbsorberMaterial = Tungsten;
  } else if ( name == "U" ||
 	      name == "Uranium" || name == "uranium" ) { 
    theAbsorberMaterial = Uranium;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Iron *  will be used." 
	   << G4endl << G4endl;  
    theAbsorberMaterial = Iron;
  }
  
  logicAbsorber->SetMaterial( theAbsorberMaterial );
  
  G4cout << G4endl
	 << " Absorber Material = " << logicAbsorber->GetMaterial()->GetName() 
	 << G4endl;
  
}


void MyDetectorConstruction::SetActiveMaterial(const G4String name) {

  if ( name == "Scintillator" || name == "scintillator" ) { 
    theActiveMaterial = Polystyrene;
  } else if ( name == "LAr" ||
 	      name == "LiquidArgon" || name == "liquidArgon" ) { 
    theActiveMaterial = LiquidArgon;
  } else if ( name == "PbWO4" ) {
    theActiveMaterial = PbWO4;
  } else if ( name == "Si" ||
 	      name == "Silicon" || name == "silicon" ) { 
    theActiveMaterial = Silicon;
  } else if ( name == "Quartz" || name == "quartz" ) { 
    theActiveMaterial = Quartz;
  } else {
    G4cout << G4endl << G4endl
	   << "WARNING: the name of the material has not been recognized!" << G4endl
	   << "     ===> the default  * Scintillator *  will be used." 
	   << G4endl << G4endl;  
    theActiveMaterial = Polystyrene;
  }
  
  logicActive->SetMaterial( theActiveMaterial );
  
  G4cout << G4endl
	 << " Active Material = " << logicActive->GetMaterial()->GetName() 
	 << G4endl;
  
}
