///////////////////////////////////////////////////////////////////////////////
// File: CCalDetectorConstruction.cc
// Description: CCalDetectorConstruction user action class to construct 
//              detector geometry
///////////////////////////////////////////////////////////////////////////////
#include "CCalDetectorConstruction.hh"

//#define debug

#ifdef debug
#include "G4Timer.hh"
#endif

#include "CCalMaterialFactory.hh"
#include "CCalRotationMatrixFactory.hh"
#include "CCalSensAssign.hh"
#include "CCalMagneticField.hh"
#include "CCalG4Hall.hh"
#include "CCalutils.hh"

#include "CCalEcalOrganization.hh"
#include "CCalHcalOrganization.hh"
#include "G4SDManager.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleRunge.hh"
#include "G4TransportationManager.hh"

CCalDetectorConstruction::CCalDetectorConstruction() {}

CCalDetectorConstruction::~CCalDetectorConstruction() {}

G4VPhysicalVolume* CCalDetectorConstruction::Construct() {

  /////////
  //Instantiate for the first time the materials and rotations
#ifdef debug
  G4cout << "Retrieving materials...." << G4endl;
#endif
  CCalMaterialFactory::getInstance("material.cms");

#ifdef debug
  G4cout << "Retrieving rotation matrices....." << G4endl;
#endif
  CCalRotationMatrixFactory::getInstance("rotation.cms");

  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------

  static G4bool fieldIsInitialized = false;
  //And finally that it was not initialized previously
  if (!fieldIsInitialized) {
    CCalMagneticField* ccalField=new CCalMagneticField("fmap.tb96");
    G4double field = ccalField->GetConstantFieldvalue();
    if (field == 0) {
      ccalField = NULL;
      G4cout << "***************************" << G4endl
	     << "*                         *" << G4endl
	     << "*  Magnetic Field is off  *" << G4endl
	     << "*                         *" << G4endl
	     << "***************************" << G4endl;
    } else {
      G4cout << "***************************" << G4endl
	     << "*                         *" << G4endl
	     << "*  Magnetic Field is on   *" << G4endl
	     << "*                         *" << G4endl
	     << "***************************" << G4endl << G4endl
	     << " Field Value " << tab << field << G4endl;
    }
    G4FieldManager* fieldMgr
      = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(ccalField);
    G4Mag_UsualEqRhs *fEquation = new G4Mag_UsualEqRhs(ccalField); 
    G4MagIntegratorStepper *pStepper = new G4SimpleRunge (fEquation);
    G4ChordFinder *pChordFinder = new G4ChordFinder(ccalField,
						    1.0e-2*mm, pStepper);
    fieldMgr->SetChordFinder(pChordFinder);
    fieldIsInitialized = true;
  }

#ifdef debug
  G4cout << tab << "CCalDetectorConstruction: Starting timer!!!" 
         << G4endl;
  G4Timer timer;
  timer.Start();
#endif

  //HCAL Test Beam 96
  CCalG4Hall*  testBeamHCal96 = new CCalG4Hall("HcalTB96");
  testBeamHCal96->constructHierarchy();
#ifdef debug
  timer.Stop();
  G4cout << tab << "CCalDetectorConstruction: Total time to "
         << "construct the geometry: " << timer << G4endl;
#endif //debug
  G4VPhysicalVolume* volume = testBeamHCal96->PhysicalVolume(0);

  //Addsenistive detector types 
  bool result;
  result = CCalSensAssign::getInstance()->addCaloSD("HadronCalorimeter",
						    new CCalHcalOrganization);
  result = CCalSensAssign::getInstance()->addCaloSD("CrystalMatrix",
						    new CCalEcalOrganization);

  //Assign the sensitive detectors
  result = CCalSensAssign::getInstance()->assign();

  //Create the stacking manager required by Calorimeter
  result = CCalSensAssign::getInstance()->stackingAction();
  
  return volume;

}
