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
// File: CCalSensitiveDetectors.hh
// Description: Container of logical volume pointers which can be sensitive
//              detectors
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalSensitiveDetectors_h
#define CCalSensitiveDetectors_h 1

#include <vector>
#include <map>
#include "G4VSensitiveDetector.hh"
#include "G4LogicalVolume.hh"

typedef std::multimap< G4String, G4LogicalVolume*, std::less<G4String> > mmslv;

class CCalSensitiveDetectors
{
public:    
  ~CCalSensitiveDetectors(){};
  std::vector<G4LogicalVolume*> getVolumes (const G4String& name,
                                            G4bool exists = false);
  void registerVolume (const G4String& name, G4LogicalVolume*);
  G4bool setSensitive(const G4String& string, G4VSensitiveDetector* sens);
  static CCalSensitiveDetectors* getInstance();
private:
  CCalSensitiveDetectors(){};
  static CCalSensitiveDetectors* theInstance;
  // logical volume container: name, G4LogicalVolume*
  mmslv theLVs;
};

#endif
