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
/// file: DamageModel.hh
/// brief: Parameters related to the DNA Damage model
#ifndef MOLECULAR_DAMAGE_MODEL
#define MOLECULAR_DAMAGE_MODEL

#include "globals.hh"
#include "DamageModelMessenger.hh"
#include "G4MolecularConfiguration.hh"
#include "G4OH.hh"
#include "G4Electron_aq.hh"
#include "G4Hydrogen.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DamageModel
{
 public:
  DamageModel();

  virtual ~DamageModel() { delete fpDamageMessenger; };

  // methods
  G4bool IsDirectStrandBreak(const G4double&) const;

  G4bool IsIndirectStrandDamage(const G4MolecularConfiguration*) const;

  G4bool IsIndirectBaseDamage(const G4MolecularConfiguration*) const;

  G4bool IsInducedStrandBreak(const G4MolecularConfiguration*) const;

  // setters
  void SetDirectDamageLower(G4double d)
  {
    fDirectDmgLower = d;
    if(fDirectDmgLower < fDirectDmgUpper)
      fDirectDmgUpper = fDirectDmgLower;
  };

  void SetDirectDamageUpper(G4double d)
  {
    fDirectDmgUpper = d;
    if(fDirectDmgLower > fDirectDmgUpper)
      fDirectDmgLower = fDirectDmgUpper;
  };

  // Indirect OH
  void SetIndirectOHStrandChance(G4double d)
  {
    CheckValidProbability("DamageStrandOH", d);
    fStrandOH = d;
  }

  void SetIndirectOHBaseChance(G4double d)
  {
    CheckValidProbability("DamageBaseOH", d);
    fBaseOH = d;
  }

  void SetInductionOHChance(G4double d)
  {
    CheckValidProbability("InductionOH", d);
    fInductionOH = d;
  };

  // Indirect Eaq
  void SetIndirectEaqStrandChance(G4double d)
  {
    CheckValidProbability("DamageStrandEaq", d);
    fStrandEaq = d;
  }

  void SetIndirectEaqBaseChance(G4double d)
  {
    CheckValidProbability("DamageBaseEaq", d);
    fBaseEaq = d;
  }

  void SetInductionEaqChance(G4double d)
  {
    CheckValidProbability("InductionEaq", d);
    fInductionEaq = d;
  };

  // Indirect H
  void SetIndirectHStrandChance(G4double d)
  {
    CheckValidProbability("DamageStrandH", d);
    fStrandH = d;
  }

  void SetIndirectHBaseChance(G4double d)
  {
    CheckValidProbability("DamageBaseH", d);
    fBaseH = d;
  }

  void SetInductionHChance(G4double d)
  {
    CheckValidProbability("InductionH", d);
    fInductionH = d;
  };

 protected:
  static void CheckValidProbability(const G4String& str, G4double p) ;

 private:
  DamageModelMessenger* fpDamageMessenger;
  G4double fDirectDmgLower = 17.5 * eV, fDirectDmgUpper = 17.5 * eV;
  G4double fInductionOH = 0, fInductionEaq = 0, fInductionH = 0;
  G4double fStrandOH = 1.0, fStrandEaq = 1.0, fStrandH = 1.0;
  G4double fBaseOH = 1.0, fBaseEaq = 1.0, fBaseH = 1.0;
  const G4MoleculeDefinition* fOH = G4OH::Definition();
  const G4MoleculeDefinition* fe_aq = G4Electron_aq::Definition();
  const G4MoleculeDefinition* fH = G4Hydrogen::Definition();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif  // MOLECULAR_DAMAGE_MODEL
