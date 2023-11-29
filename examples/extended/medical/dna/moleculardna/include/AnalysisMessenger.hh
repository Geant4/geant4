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
#ifndef MOLECULAR_DETECTOR_MESSENGER_HH
#define MOLECULAR_DETECTOR_MESSENGER_HH

#include "globals.hh"
#include "G4UImessenger.hh"
#include <memory>
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

class AnalysisManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class AnalysisMessenger : public G4UImessenger
{
 public:
  explicit AnalysisMessenger(AnalysisManager*);

  ~AnalysisMessenger() override = default;

  void SetNewValue(G4UIcommand*, G4String) override;

 protected:
 private:
  AnalysisManager* fpAnalysisManager;

  // Related to geometry
  std::unique_ptr<G4UIdirectory> fpAnalysisDirectory;
  std::unique_ptr<G4UIcmdWithABool> fpSaveStrands;
  std::unique_ptr<G4UIcmdWithAString> fpStrandDirectory;
  std::unique_ptr<G4UIcmdWithAnInteger> fpFragmentLength;
  std::unique_ptr<G4UIcmdWithAnInteger> fpSaveSingleChain;
  std::unique_ptr<G4UIcmdWithAnInteger> fpDSBDistance;
  std::unique_ptr<G4UIcmdWithoutParameter> fpTestClassifier;
  std::unique_ptr<G4UIcmdWithAString> fpFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_DETECTOR_MESSENGER_HH
