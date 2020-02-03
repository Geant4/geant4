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
/// \file field/field03/include/F03CalorimeterSD.hh
/// \brief Definition of the F03CalorimeterSD class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F03CalorimeterSD_h
#define F03CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "F03CalorHit.hh"

class F03DetectorConstruction;
class G4HCofThisEvent;
class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F03CalorimeterSD : public G4VSensitiveDetector
{
  public:

      F03CalorimeterSD(G4String, F03DetectorConstruction* );
      virtual ~F03CalorimeterSD();

      virtual void Initialize(G4HCofThisEvent*);
      virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      virtual void EndOfEvent(G4HCofThisEvent*);

  private:

      F03CalorHitsCollection*  fCalCollection;
      F03DetectorConstruction* fDetector;
      G4int*                   fHitID;
};

#endif
