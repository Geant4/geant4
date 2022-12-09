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

#ifndef eRositaTrackerSD_h
#define eRositaTrackerSD_h 1

#include "eRositaTrackerHit.hh"

#include "G4VSensitiveDetector.hh"

class G4HCofThisEvent;
class G4Step;

class eRositaTrackerSD : public G4VSensitiveDetector {
public:    
    explicit eRositaTrackerSD(G4String name);
    
    ~eRositaTrackerSD() override;

    void EndOfEvent(G4HCofThisEvent* collection) override;
    
    void Initialize(G4HCofThisEvent* collection) override;
    
    auto ProcessHits(G4Step* step, G4TouchableHistory* history) -> G4bool override;

private:
    eRositaTrackerHitsCollection* trackerCollection;
};
#endif
