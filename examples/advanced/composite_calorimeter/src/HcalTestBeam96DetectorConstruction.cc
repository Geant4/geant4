///////////////////////////////////////////////////////////////////////////////
// File: HcalTestBeam96DetectorConstruction.cc
// Author: S.Banerjee
// Modifications: 23/08/00 
///////////////////////////////////////////////////////////////////////////////
#include "HcalTestBeam96DetectorConstruction.hh"

//#define debug

#ifdef debug
#include "G4Timer.hh"
#endif

#include "CCalMaterialFactory.hh"
#include "CCalRotationMatrixFactory.hh"
#include "CCalSensAssign.hh"
#include "TestBeamMagneticField.hh"
#include "G4HcalTB96.hh"
#include "utils.hh"

#include "G4CaloSD.hh"
#include "CCalEcalOrganization.hh"
#include "CCalHcalOrganization.hh"
#include "G4SDManager.hh"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleRunge.hh"
#include "G4TransportationManager.hh"

HcalTestBeam96DetectorConstruction::HcalTestBeam96DetectorConstruction() {}

HcalTestBeam96DetectorConstruction::~HcalTestBeam96DetectorConstruction() {}

G4VPhysicalVolume* HcalTestBeam96DetectorConstruction::Construct() {

  /////////
  //Instantiate for the first time the materials and rotations
#ifdef debug
  cout << "Retrieving materials...." << endl;
#endif
  CCalMaterialFactory::getInstance("material.cms");

#ifdef debug
  cout << "Retrieving rotation matrices....." << endl;
#endif
  CCalRotationMatrixFactory::getInstance("rotation.cms");

  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------

  static G4bool fieldIsInitialized = false;
  //And finally that it was not initialized previously
  if (!fieldIsInitialized) {
    TestBeamMagneticField* ccalField=new TestBeamMagneticField("fmap.tb96");
    G4double field = ccalField->GetConstantFieldvalue();
    if (field == 0) {
      ccalField = NULL;
      cout << "***************************" << endl
           << "*                         *" << endl
           << "*  Magnetic Field is off  *" << endl
           << "*                         *" << endl
	   << "***************************" << endl;
    } else {
      cout << "***************************" << endl
           << "*                         *" << endl
           << "*  Magnetic Field is on   *" << endl
           << "*                         *" << endl
	   << "***************************" << endl << endl
	   << " Field Value " << tab << field << endl;
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
  G4cout << tab << "HcalTestBeam96DetectorConstruction: Starting timer!!!" 
         << endl;
  G4Timer timer;
  timer.Start();
#endif

  //HCAL Test Beam 96
  G4HcalTB96*  testBeamHCal96 = new G4HcalTB96("HcalTB96");
  testBeamHCal96->constructHierarchy();
#ifdef debug
  timer.Stop();
  G4cout << tab << "HcalTestBeam96DetectorConstruction: Total time to "
         << "construct the geometry: " << timer << endl;
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
