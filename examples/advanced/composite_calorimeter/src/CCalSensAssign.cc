///////////////////////////////////////////////////////////////////////////////
// File: CCalSensAssign.cc
// Description: CCalSenAssign creates and assigns the sensitive detetctors 
///////////////////////////////////////////////////////////////////////////////
#include "CCalSensAssign.hh"
#include "CCalStackingAction.hh"
#include "CCalSensitiveDetectors.hh"
#include "CCaloSD.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"

//#define debug

CCalSensAssign* CCalSensAssign::theInstance = 0;


CCalSensAssign* CCalSensAssign::getInstance() {
 if (!theInstance) theInstance = new CCalSensAssign;
 return theInstance;
}


CCalSensAssign::CCalSensAssign() {}


bool CCalSensAssign::assign() {
  bool result = false;

  CCalSensitiveDetectors* sensDets = CCalSensitiveDetectors::getInstance();
  for (G4std::map<G4String,G4VSensitiveDetector*>::const_iterator sens_it = sens_.begin();
       sens_it!=sens_.end(); ++sens_it) {
    G4String name = sens_it->first;
    G4VSensitiveDetector* sens = sens_it->second;
    if (sensDets->setSensitive(name, sens)) {
      G4SDManager::GetSDMpointer()->AddNewDetector(sens);
#ifdef debug
      G4cout << "Add " << sens->GetName() 
	   << " to the list of Sensitive detetctors" << G4endl;
#endif
    }
  }

  return result;
}

bool CCalSensAssign::stackingAction() {

  bool result = false;
  //Create the stacking manager required by Calorimeter
  if (G4RunManager::GetRunManager()->GetUserStackingAction() == 0) {
    G4cout << "***CCalSensAssign creating a CCalStackingAction ***" << G4endl;
    G4RunManager::GetRunManager()->SetUserAction(new CCalStackingAction);
    result = true;
  } else {
    G4cout << "***CCalSens: a StackingAction already exists. "
	 << "Maybe not the one CCaloSD needs?" << G4endl;
  }  
  return result;
}


bool CCalSensAssign::addCaloSD(G4String name, 
			       CCalVOrganization* numberingScheme) {
  sens_[name]         = new CCaloSD(name, numberingScheme);
  return true;

}
