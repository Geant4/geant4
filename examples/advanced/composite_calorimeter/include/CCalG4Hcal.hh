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
// File: CCalG4Hcal.hh
// Description: Euipped to construct the G4 geometry of the hadron calorimeter
//              in the 96 test beam run
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Hcal_h
#define CCalG4Hcal_h 1

#include "CCalHcal.hh"
#include "CCalG4Able.hh"
#include <vector>

typedef G4LogicalVolume* ptrG4Log;

class CCalG4Hcal: public CCalHcal, public CCalG4Able {
public:
  //Constructor and Destructor
  CCalG4Hcal(const G4String &name);
  virtual ~CCalG4Hcal();

protected:
  //This methods actually constructs the volume.
  virtual G4VPhysicalVolume* constructIn(G4VPhysicalVolume*);
  virtual void constructDaughters();

  //Construct layer of scintillator or absorber
  G4LogicalVolume*   constructScintillatorLayer (G4int);
  G4LogicalVolume*   constructAbsorberLayer (G4int);

  //Constructs the sensitive detectors and associates them to the corresponding
  //logical volumes
  virtual void constructSensitive() ;

private:
  //Private data members
  ptrG4Log*        sclLog;
  ptrG4Log*        absLog;

  // Logical volumes for sensitive detectors
  std::vector<ptrG4Log> allSensitiveLogs;
};

#endif
