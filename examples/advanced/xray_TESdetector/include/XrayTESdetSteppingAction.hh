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
/// \file XrayTESdetSteppingAction.hh
/// \brief Definition of the SteppingAction class
//
//
// Authors: P.Dondero (paolo.dondero@cern.ch), R.Stanzani (ronny.stanzani@cern.ch)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#ifndef XrayTESdetSteppingAction_h
#define XrayTESdetSteppingAction_h 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"
#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class XrayTESdetDetectorConstruction;
class XrayTESdetRunAction;
class XrayTESdetEventAction;

class XrayTESdetSteppingAction : public G4UserSteppingAction
{
  public:

    explicit XrayTESdetSteppingAction() = default;
    ~XrayTESdetSteppingAction() override = default;
    void UserSteppingAction(const G4Step*) override;

  private:

    G4int fPrev_eventID = 0;
    G4int fPrev_trackID = 0;
    std::map<G4int, G4double> fInit_energy;
    std::map<G4int, G4String> fCreator_proc;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif