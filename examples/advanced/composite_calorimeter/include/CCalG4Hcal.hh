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
// File: CCalG4Hcal.hh
// Description: Euipped to construct the G4 geometry of the hadron calorimeter
//              in the 96 test beam run
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalG4Hcal_h
#define CCalG4Hcal_h 1

#include "CCalHcal.hh"
#include "CCalG4Able.hh"
#include "g4std/vector"

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
  G4std::vector<ptrG4Log> allSensitiveLogs;
};

#endif
