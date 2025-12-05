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
/// \file GB03EventAction.hh
/// \brief Definition of the GB03EventAction class

#ifndef GB03EventAction_h
#define GB03EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
template<typename T>
class G4THitsMap;

class GB03EventAction : public G4UserEventAction
{
  public:
    GB03EventAction() = default;
    ~GB03EventAction() override = default;

    void BeginOfEventAction(const G4Event* event) override;
    void EndOfEventAction(const G4Event* event) override;

  private:
    // methods
    G4THitsMap<G4double>* GetHitsCollection(G4int hcID, const G4Event* event) const;
    G4double GetSum(G4THitsMap<G4double>* hitsMap) const;
    void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength, G4double gapEdep,
                              G4double gapTrackLength) const;

    // data members
    G4int fAbsoEdepHCID = -1;
    G4int fGapEdepHCID = -1;
    G4int fAbsoTrackLengthHCID = -1;
    G4int fGapTrackLengthHCID = -1;
};

#endif
