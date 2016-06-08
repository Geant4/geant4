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
//
///////////////////////////////////////////////////////////////////////////////
// File: CCalSensitiveDetectors.hh
// Description: Container of logical volume pointers which can be sensitive
//              detectors
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensitiveDetectors_h
#define CCalSensitiveDetectors_h 1

#include "g4std/vector"
#include "g4std/map"
#include "G4VSensitiveDetector.hh"
#include "G4LogicalVolume.hh"

typedef G4std::multimap< G4String, G4LogicalVolume*, G4std::less<G4String> > mmslv;

class CCalSensitiveDetectors {

public:    
  ~CCalSensitiveDetectors(){};
  G4std::vector<G4LogicalVolume*> getVolumes (const G4String& name, bool exists = 0);
  void registerVolume (const G4String& name, G4LogicalVolume*);
  bool setSensitive(const G4String& string, G4VSensitiveDetector* sens);
  static CCalSensitiveDetectors* getInstance();
private:
  CCalSensitiveDetectors(){};
  static CCalSensitiveDetectors* theInstance;
  // logical volume container: name, G4LogicalVolume*
  mmslv theLVs;

};

#endif
