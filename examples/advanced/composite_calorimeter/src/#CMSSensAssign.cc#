#include "CMSSensAssign.hh"
#include "CMSStackingAction.hh"
#include "CMSSensitiveDetectors.hh"
#include "G4CaloSD.hh"

#include "G4RunManager.hh"
#include "G4SDManager.hh"

//#define debug

CMSSensAssign* CMSSensAssign::theInstance = 0;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CMSSensAssign* CMSSensAssign::getInstance() {

 if (!theInstance) theInstance = new CMSSensAssign;
 return theInstance;
}

CMSSensAssign::CMSSensAssign() {}


bool CMSSensAssign::assign() {
  bool result = false;

  CMSSensitiveDetectors* sensDets = CMSSensitiveDetectors::getInstance();
  for (map<G4String,G4VSensitiveDetector*>::const_iterator sens_it = sens_.begin();
       sens_it!=sens_.end(); ++sens_it) {
    G4String name = sens_it->first;
    G4VSensitiveDetector* sens = sens_it->second;
    //    cout << " CMSSensAssign " << name << " Sens " << sens << endl;
    //    vector<G4LogicalVolume*> lvs = sensDets->getVolumes (name, 1);
    //    if (lvs.size()>0) {
    if (sensDets->setSensitive(name, sens)) {
      G4SDManager::GetSDMpointer()->AddNewDetector(sens);
      //      for (vector<G4LogicalVolume*>::iterator iter=lvs.begin(); iter<lvs.end();
      //	   iter++) {
      //	(*iter)->SetSensitiveDetector(sens);
#ifdef debug
      //      cout << " Associate SD " << name << " to " << (*iter)->GetName() 
      //     << endl;
      cout << "Add " << sens->GetName() 
	   << " to the list of Sensitive detetctors" << endl;
#endif
      //      }
    }
  }

  return result;
}

bool CMSSensAssign::stackingAction() {

  bool result = false;
  //Create the stacking manager required by Calorimeter
  if (G4RunManager::GetRunManager()->GetUserStackingAction() == 0) {
    cout << "***CMSSensAssign creating a CMSStackingAction ***" << endl;
    G4RunManager::GetRunManager()->SetUserAction(new CMSStackingAction);
    result = true;
  } else {
    cout << "***CMSSens: a StackingAction already exists. "
	 << "Maybe not the one G4CaloSD needs?" << endl;
  }  
  return result;
}

bool CMSSensAssign::addCaloSD(G4String name, 
			      VDetectorOrganization* numberingScheme) {
  sens_[name]         = new G4CaloSD(name, numberingScheme);
  return true;

}
