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
///////////////////////////////////////////////////////////////////////////////

#include "CCalSensitiveDetectors.hh"

//#define debug

CCalSensitiveDetectors* CCalSensitiveDetectors::theInstance = 0;


CCalSensitiveDetectors* CCalSensitiveDetectors::getInstance() {
 if (!theInstance) theInstance = new CCalSensitiveDetectors;
 return theInstance;
}


void CCalSensitiveDetectors::registerVolume (const G4String& string, 
					     G4LogicalVolume* logv) {

  theLVs.insert(mmslv::value_type(string, logv));
#ifdef debug
  G4cout << "CCalSensitiveDetectors : Register " << logv->GetName() 
       << " in category " << string << G4endl;
#endif
}

G4std::vector<G4LogicalVolume*> CCalSensitiveDetectors::getVolumes (const G4String& string, 
							     bool exist) {

  mmslv::const_iterator mmscite;
  G4std::pair<mmslv::iterator, mmslv::iterator> mmsdi;
  mmsdi = theLVs.equal_range(string);
  G4std::vector<G4LogicalVolume*> lvs;
  for (mmscite = mmsdi.first; mmscite != mmsdi.second; mmscite++ ) {
    lvs.push_back(const_cast<G4LogicalVolume*>((*mmscite).second));
  }

  if (exist) G4cout << "CCalSensitiveDetector : " << lvs.size() 
		  << " detectors for " << string << G4endl;
  
  return lvs;
}


bool CCalSensitiveDetectors::setSensitive(const G4String& string, 
					  G4VSensitiveDetector* sens) {

  bool result=false;
  mmslv::const_iterator mmscite;
  G4std::pair<mmslv::iterator, mmslv::iterator> mmsdi;
  mmsdi = theLVs.equal_range(string);
  for (mmscite = mmsdi.first; mmscite != mmsdi.second; mmscite++ ) {
    G4LogicalVolume* lv = const_cast<G4LogicalVolume*>((*mmscite).second);
    lv ->SetSensitiveDetector(sens);
    result = true;
#ifdef debug
    G4cout << " Associate SD " << sens->GetName() << " to " << lv->GetName() 
	 << G4endl;
#endif
  }
  return result;

}
