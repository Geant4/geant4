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
/// \file SAXSSensitiveDetector.hh
/// \brief Definition of the SAXSSensitiveDetector class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SAXSSensitiveDetector_h
#define SAXSSensitiveDetector_h 1

#include "SAXSSensitiveDetectorHit.hh"
#include "SAXSSensitiveDetectorMessenger.hh"
#include "G4VSensitiveDetector.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/// Sensitive detector.
/// It stores the features of the impinging particles.

class SAXSSensitiveDetector : public G4VSensitiveDetector
{
public:
    SAXSSensitiveDetector(G4String);
    virtual ~SAXSSensitiveDetector();
    
    virtual void Initialize(G4HCofThisEvent*);
    virtual G4bool ProcessHits(G4Step*,G4TouchableHistory*);
    virtual void EndOfEvent(G4HCofThisEvent*);
    
    void SetVarStopAndKill(G4bool bVar) {fVarStopAndKill = bVar;}
    
private:
    SensitiveDetectorHitsCollection* fHitsCollection;
    G4int fHCID;  
    
    SAXSSensitiveDetectorMessenger* fSDMessenger;
    
    G4bool fVarStopAndKill;    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

