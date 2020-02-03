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
///////////////////////////////////////////////////////////////////////////////
// File: CCalSensitiveDetectors.cc
// Description: Provides a map of logicalvolume pointers which can be sensitive
///////////////////////////////////////////////////////////////////////////////

#include "CCalSensitiveDetectors.hh"

//#define debug

CCalSensitiveDetectors* CCalSensitiveDetectors::theInstance = 0;


CCalSensitiveDetectors* CCalSensitiveDetectors::getInstance()
{
 if (!theInstance) theInstance = new CCalSensitiveDetectors;
 return theInstance;
}


void CCalSensitiveDetectors::registerVolume (const G4String& string, 
                                             G4LogicalVolume* logv)
{
  theLVs.insert(mmslv::value_type(string, logv));
#ifdef debug
  G4cout << "CCalSensitiveDetectors : Register " << logv->GetName() 
       << " in category " << string << G4endl;
#endif
}

std::vector<G4LogicalVolume*>
CCalSensitiveDetectors::getVolumes (const G4String& string, G4bool exist)
{
  mmslv::const_iterator mmscite;
  std::pair<mmslv::iterator, mmslv::iterator> mmsdi;
  mmsdi = theLVs.equal_range(string);
  std::vector<G4LogicalVolume*> lvs;
  for (mmscite = mmsdi.first; mmscite != mmsdi.second; mmscite++ ) {
    lvs.push_back(const_cast<G4LogicalVolume*>((*mmscite).second));
  }

  if (exist) G4cout << "CCalSensitiveDetector : " << lvs.size() 
                  << " detectors for " << string << G4endl;
  
  return lvs;
}


G4bool CCalSensitiveDetectors::setSensitive(const G4String& string, 
                                          G4VSensitiveDetector* sens)
{
  G4bool result=false;
  mmslv::const_iterator mmscite;
  std::pair<mmslv::iterator, mmslv::iterator> mmsdi;
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
