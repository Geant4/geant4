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
//
/// \file EventAction.hh
/// \brief Definition of the EventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

namespace B4a
{

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    ~EventAction() override;

    void  BeginOfEventAction(const G4Event* event) override;
    void    EndOfEventAction(const G4Event* event) override;

    void AddAbs(G4double de, G4double dl);
    void AddGap(G4double de, G4double dl);

  private:
    G4double  fEnergyAbs = 0.;
    G4double  fEnergyGap = 0.;
    G4double  fTrackLAbs = 0.;
    G4double  fTrackLGap = 0.;
};

// inline functions

inline void EventAction::AddAbs(G4double de, G4double dl) {
  fEnergyAbs += de;
  fTrackLAbs += dl;
}

inline void EventAction::AddGap(G4double de, G4double dl) {
  fEnergyGap += de;
  fTrackLGap += dl;
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


