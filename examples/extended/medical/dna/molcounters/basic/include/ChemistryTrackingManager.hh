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
/// \file ChemistryTrackingManager.hh
/// \brief Definition of the ChemistryTrackingManager class

// The `molcounters` example(s) are provided as part of Geant4-DNA
// and any report or published result obtained using it shall cite
// the respective Geant4-DNA collaboration publications.
//
// Reports or results obtained using the spatially-aware `MoleculeCounter`
// provided in this example, shall further cite:
//
// Velten & Tomé, Radiation Physics and Chemistry, 2023 (10.1016/j.radphyschem.2023.111194)
//
//
// Author: Christian Velten (2025)
//

#ifndef ChemistryTrackingManager_hh
#define ChemistryTrackingManager_hh 1

#include "G4ITTrackingInteractivity.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ChemistryTrackingManager : public G4ITTrackingInteractivity
{
  public:
    ChemistryTrackingManager() = default;
    ~ChemistryTrackingManager() override;

    void AppendStep(G4Track*, G4Step*) override;
    void Finalize() override;

    void SetUserAction(G4UserSteppingAction*);
    const G4UserSteppingAction* GetUserSteppingAction() const;

  private:
    G4UserSteppingAction* fUserSteppingAction = nullptr;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline void ChemistryTrackingManager::SetUserAction(G4UserSteppingAction* steppingAction)
{
  fUserSteppingAction = steppingAction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
inline const G4UserSteppingAction* ChemistryTrackingManager::GetUserSteppingAction() const
{
  return fUserSteppingAction;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#endif
