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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronSensitiveFoil.hh
//
#ifndef STCyclotronSensitiveFoil_h
#define STCyclotronSensitiveFoil_h 1

#include "G4VSensitiveDetector.hh"
#include "STCyclotronDetectorConstruction.hh"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class STCyclotronRun;

class STCyclotronSensitiveFoil : public G4VSensitiveDetector
{

public:

  STCyclotronSensitiveFoil(const G4String& name,
			   STCyclotronDetectorConstruction* det);
  virtual ~STCyclotronSensitiveFoil();
  
  //methods from base class
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
 
private:
  STCyclotronDetectorConstruction* fDet;
  STCyclotronRun* fRun;
  G4int fTempTrack;
  G4int fTempTrack1;
  G4ThreeVector fTempVector;
  G4double fTempEnergy;

};

#endif
