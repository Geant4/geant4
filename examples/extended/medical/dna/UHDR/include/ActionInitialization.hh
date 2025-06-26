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
/// \file ActionInitialization.hh
/// \brief Definition of the ActionInitialization class

#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4GenericMessenger.hh"
#include "G4String.hh"
#include "G4VUserActionInitialization.hh"

#include <memory>



class ActionInitialization : public G4VUserActionInitialization
{
  public:
    explicit ActionInitialization();

    ~ActionInitialization() override = default;

    void BuildForMaster() const override;

    void Build() const override;

  private:
    void SetPulseStructureHistoInput(G4String);  // L.T.Anh: set histo file
    void SetPulseStructureInput(G4String);
    void SetPulsePeriod(G4double tp) { fPulsePeriod = tp; }
    void SetNumberOfPulse(G4int npulse) { fNumberOfPulse = npulse; }
    void DefineCommands();  // Le Tuan Anh: add new commands
   
    G4String fPulseStructure = "";
    // Le Tuan Anh: To store pulse structure filename set by user in macro:
    std::unique_ptr<G4GenericMessenger> fMessenger;  // Le Tuan Anh: command control
    G4bool fActivePulse = false;  // Le Tuan Anh: flag for invoking pulse mode
    G4bool fUseHistoInput = false;
    G4bool fUseInterPulse = false;  // Le Tuan Anh: unify with fActivePulse later
    G4double fPulsePeriod = 0;
    G4int fNumberOfPulse = 1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
