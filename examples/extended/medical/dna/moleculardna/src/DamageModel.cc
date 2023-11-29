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
/// file: DamageModel.cc
/// brief: Parameters related to the DNA Damage model

#include "DamageModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DamageModel::DamageModel()
{
  fpDamageMessenger = new DamageModelMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DamageModel::IsDirectStrandBreak(const G4double& d) const
{
  G4bool bool_val = false;
  if(d < fDirectDmgLower)
  {
    bool_val = false;
  }
  else if(d >= fDirectDmgUpper)
  {
    bool_val = true;
  }
  else 
  {
    bool_val = (G4UniformRand() < ((d - fDirectDmgLower) / (fDirectDmgUpper - fDirectDmgLower)));
  }
  return bool_val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DamageModel::IsIndirectBaseDamage(const G4MolecularConfiguration* mol) const
{
  G4bool bool_val = false;
  if(mol->GetDefinition() == fOH)
  {
    bool_val = (G4UniformRand() < fBaseOH); 
  }
  else if(mol->GetDefinition() == fH)
  {
    bool_val = (G4UniformRand() < fBaseH);
  }
  else if(mol->GetDefinition() == fe_aq)
  {
    bool_val = (G4UniformRand() < fBaseEaq);
  }
  else
  {
    bool_val = false;
  }
  return bool_val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DamageModel::IsIndirectStrandDamage(const G4MolecularConfiguration* mol) const
{
  G4bool bool_val = false;
  if(mol->GetDefinition() == fOH)
  {
    bool_val = (G4UniformRand() < fStrandOH);
  }
  else if(mol->GetDefinition() == G4Hydrogen::Definition())
  {
    bool_val = (G4UniformRand() < fStrandH);
  }
  else if(mol->GetDefinition() == G4Electron_aq::Definition())
  {
    bool_val = (G4UniformRand() < fStrandEaq);
  }
  else
  {
    bool_val = false;
  }
  return bool_val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DamageModel::IsInducedStrandBreak(
    const G4MolecularConfiguration* mol) const
{
  G4bool bool_val = false;
  if(mol->GetDefinition() == fOH)
  {
    bool_val = (G4UniformRand() < fInductionOH);
  }
  else if(mol->GetDefinition() == G4Hydrogen::Definition())
  {
    bool_val = (G4UniformRand() < fInductionH);
  }
  else if(mol->GetDefinition() == G4Electron_aq::Definition())
  {
    bool_val = (G4UniformRand() < fInductionEaq);
  }
  else
  {
    bool_val = false;
  }
  return bool_val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DamageModel::CheckValidProbability(const G4String& str, G4double p)
{
  if((p < 0) || (p > 1))
  {
    G4ExceptionDescription errmsg;
    errmsg << "Invalid probability for " << str;
    G4Exception("DamageModel", "ERR_INVALID_PROB", FatalException, errmsg);
  }
}