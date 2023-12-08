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
/// \file PhysSteppingAction.hh
/// \brief Definition of the PhysSteppingAction class

#ifndef PhysSteppingAction_h
#define PhysSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "PhysEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class PhysSteppingAction : public G4UserSteppingAction
{
public:
    PhysSteppingAction(PhysEventAction* pEvent);
    ~PhysSteppingAction() override = default;

    void UserSteppingAction(const G4Step*step) override; 

    G4double SetupVolumeFlag(const G4String &volumeName);
    G4double SetupProcessFlag(const G4String &processName);
    G4double SetupParticleFlag(const G4String &particleName);

private:
    PhysEventAction* fEventAction{nullptr};
    G4double fFlagProcess{0.};
    G4double fFlagVolume{0.};
    G4double fFlagParticle{0.};
    G4double fFlagParentID{0.};
    G4double fLastMetVoxelCopyNumber{-1.};

    void SetupFlags(const G4Step *step);
    void SetupVoxelCopyNumber(const G4Step* step);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
