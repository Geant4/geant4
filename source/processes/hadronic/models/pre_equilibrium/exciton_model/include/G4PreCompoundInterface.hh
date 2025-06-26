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
// GEANT4 Class header file
//
// File name:  G4PreCompoundInterface
//
// Author:  V.Ivantchenko, 20 January 2025
//
// Class Description:
// Model implementation for pre-equilibrium decay model.
// It is an alternative to the default model.
//

#ifndef G4PreCompoundInterface_h
#define G4PreCompoundInterface_h 1

#include "G4VPreCompoundModel.hh"
#include "G4Fragment.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ExcitationHandler.hh"

class G4PreCompoundEmissionInt;
class G4VPreCompoundTransitions;
class G4NuclearLevelData;
class G4ParticleDefinition;

class G4PreCompoundInterface : public G4VPreCompoundModel
{ 
public:

  G4PreCompoundInterface();

  ~G4PreCompoundInterface() override;

  G4ReactionProductVector* DeExcite(G4Fragment& aFragment) override;

  void BuildPhysicsTable(const G4ParticleDefinition&) override;
  void InitialiseModel() override;

  void ModelDescription(std::ostream& outFile) const override;
  void DeExciteModelDescription(std::ostream&) const override;
  
  G4PreCompoundInterface(const G4PreCompoundInterface &) = delete;
  const G4PreCompoundInterface& operator=
  (const G4PreCompoundInterface &right) = delete;
  G4bool operator==(const G4PreCompoundInterface &right) const = delete;
  G4bool operator!=(const G4PreCompoundInterface &right) const = delete;

private:  

  void BreakUpFragment(G4Fragment&, G4ReactionProductVector*);

  inline 
  void PerformEquilibriumEmission(const G4Fragment&, 
				  G4ReactionProductVector*) const;

  G4PreCompoundEmissionInt* theEmission{nullptr};
  G4VPreCompoundTransitions* theTransition{nullptr};
  G4NuclearLevelData* fNuclData{nullptr};

  G4double fLowLimitExc{0.0};
  G4double fHighLimitExc{DBL_MAX};

  G4bool isActive{true};
  G4bool isInitialised{false};

  G4int minZ{9};
  G4int minA{17};
  G4int fVerbose{1};
};

inline void G4PreCompoundInterface::PerformEquilibriumEmission(
            const G4Fragment& aFragment,
            G4ReactionProductVector* result) const 
{
  auto deexResult = GetExcitationHandler()->BreakItUp(aFragment);
  for (auto & frag : *deexResult) { result->push_back(std::move(frag)); }
  delete deexResult;
}

#endif

