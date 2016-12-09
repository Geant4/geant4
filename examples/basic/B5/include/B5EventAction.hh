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
// $Id: B5EventAction.hh 101036 2016-11-04 09:00:23Z gcosmo $
//
/// \file B5EventAction.hh
/// \brief Definition of the B5EventAction class

#ifndef B5EventAction_h
#define B5EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

/// Event action

class B5EventAction : public G4UserEventAction
{
public:
    B5EventAction();
    virtual ~B5EventAction();
    
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    std::vector<G4double>& GetEmCalEdep() { return fEmCalEdep; }
    std::vector<G4double>& GetHadCalEdep() { return fEmCalEdep; }
    
private:
    G4int fHodHC1ID;
    G4int fHodHC2ID;
    G4int fDriftHC1ID;
    G4int fDriftHC2ID;
    G4int fEmCalHCID;
    G4int fHadCalHCID;
    std::vector<G4double> fEmCalEdep;
    std::vector<G4double> fHadCalEdep;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
