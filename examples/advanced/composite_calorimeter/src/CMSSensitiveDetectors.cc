///////////////////////////////////////////////////////////////////////////////
// File: CMSSensitiveDetectors.cc
// Date: 15/05/02
// Modifications: 
///////////////////////////////////////////////////////////////////////////////

#include "CMSSensitiveDetectors.hh"

//#define debug

CMSSensitiveDetectors* CMSSensitiveDetectors::theInstance = 0;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CMSSensitiveDetectors* CMSSensitiveDetectors::getInstance() {

 if (!theInstance) theInstance = new CMSSensitiveDetectors;
 return theInstance;
}


void CMSSensitiveDetectors::registerVolume (const string& string, G4LogicalVolume* logv) {

  theLVs.insert(mmslv::value_type(string, logv));
#ifdef debug
  cout << "CMSSensitiveDetectors : Register " << logv->GetName() 
       << " in category " << string << endl;
#endif
}

vector<G4LogicalVolume*> CMSSensitiveDetectors::getVolumes (const string& string, bool exist) {

  mmslv::const_iterator mmscite;
  pair<mmslv::iterator, mmslv::iterator> mmsdi;
  mmsdi = theLVs.equal_range(string);
  vector<G4LogicalVolume*> lvs;
  for (mmscite = mmsdi.first; mmscite != mmsdi.second; mmscite++ ) {
    lvs.push_back(const_cast<G4LogicalVolume*>((*mmscite).second));
  }

  if (exist) cout << "CMSSensitiveDetector : " << lvs.size() 
		  << " detectors for " << string << endl;
  
  return lvs;
}


bool CMSSensitiveDetectors::setSensitive(const string& string, G4VSensitiveDetector* sens) {

  bool result=false;
  mmslv::const_iterator mmscite;
  pair<mmslv::iterator, mmslv::iterator> mmsdi;
  mmsdi = theLVs.equal_range(string);
  for (mmscite = mmsdi.first; mmscite != mmsdi.second; mmscite++ ) {
    G4LogicalVolume* lv = const_cast<G4LogicalVolume*>((*mmscite).second);
    lv ->SetSensitiveDetector(sens);
    result = true;
#ifdef debug
    cout << " Associate SD " << sens->GetName() << " to " << lv->GetName() 
	 << endl;
#endif
  }
  return result;

}
