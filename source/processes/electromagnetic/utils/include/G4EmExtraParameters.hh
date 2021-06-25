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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4EmExtraParameters
//
// Author:        Vladimir Ivanchenko
//                  
// Creation date: 06.05.2019
//
// Class Description:
//
// An internal utility class, responsable for keeping parameters
// for EM processes and models.
//
// It is initialized by the master thread but can be updated 
// at any moment via G4EmParameters interface. It is not assumed
// to be used for a direct initialisation
//
// -------------------------------------------------------------------
//

#ifndef G4EmExtraParameters_h
#define G4EmExtraParameters_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4EmExtraParametersMessenger;
class G4VEnergyLossProcess;
class G4VEmProcess;
class G4ParticleDefinition;
class G4VAtomDeexcitation;

class G4EmExtraParameters
{
public:

  explicit G4EmExtraParameters();

  ~G4EmExtraParameters();

  void Initialise();

  G4bool GetDirectionalSplitting();
  void SetDirectionalSplitting(G4bool v);

  G4bool QuantumEntanglement();
  void SetQuantumEntanglement(G4bool v);

  void SetDirectionalSplittingRadius(G4double r);
  G4double GetDirectionalSplittingRadius();

  void SetDirectionalSplittingTarget(const G4ThreeVector& v);
  G4ThreeVector GetDirectionalSplittingTarget() const;

  void SetStepFunction(G4double v1, G4double v2);
  G4double GetStepFunctionP1() const;
  G4double GetStepFunctionP2() const;
 
  void SetStepFunctionMuHad(G4double v1, G4double v2);
  G4double GetStepFunctionMuHadP1() const;
  G4double GetStepFunctionMuHadP2() const;

  void SetStepFunctionLightIons(G4double v1, G4double v2);
  G4double GetStepFunctionLightIonsP1() const;
  G4double GetStepFunctionLightIonsP2() const;

  void SetStepFunctionIons(G4double v1, G4double v2);
  G4double GetStepFunctionIonsP1() const;
  G4double GetStepFunctionIonsP2() const;

  void FillStepFunction(const G4ParticleDefinition*, G4VEnergyLossProcess*) const;

  // parameters per region or per process 
  void AddPAIModel(const G4String& particle,
                   const G4String& region,
                   const G4String& type);
  const std::vector<G4String>& ParticlesPAI() const;
  const std::vector<G4String>& RegionsPAI() const;
  const std::vector<G4String>& TypesPAI() const;

  void AddPhysics(const G4String& region, const G4String& type);
  const std::vector<G4String>& RegionsPhysics() const;
  const std::vector<G4String>& TypesPhysics() const;

  void SetSubCutRegion(const G4String& region);

  void SetProcessBiasingFactor(const G4String& procname, 
                               G4double val, G4bool wflag);

  void ActivateForcedInteraction(const G4String& procname, 
                                 const G4String& region,
                                 G4double length, 
                                 G4bool wflag);

  void ActivateSecondaryBiasing(const G4String& name,
				const G4String& region, 
				G4double factor,
				G4double energyLimit);

  // initialisation methods
  void DefineRegParamForLoss(G4VEnergyLossProcess*) const;
  void DefineRegParamForEM(G4VEmProcess*) const;

  G4EmExtraParameters(G4EmExtraParameters &) = delete;
  G4EmExtraParameters & operator=
  (const G4EmExtraParameters &right) = delete;  

private:

  G4String CheckRegion(const G4String&) const;

  void PrintWarning(G4ExceptionDescription& ed) const;

  G4EmExtraParametersMessenger* theMessenger;

  G4bool directionalSplitting;
  G4bool quantumEntanglement;

  G4double dRoverRange;
  G4double finalRange;
  G4double dRoverRangeMuHad;
  G4double finalRangeMuHad;
  G4double dRoverRangeLIons;
  G4double finalRangeLIons;
  G4double dRoverRangeIons;
  G4double finalRangeIons;

  G4double directionalSplittingRadius;
  G4ThreeVector directionalSplittingTarget;

  std::vector<G4String>  m_particlesPAI;
  std::vector<G4String>  m_regnamesPAI;
  std::vector<G4String>  m_typesPAI;

  std::vector<G4String>  m_regnamesPhys;
  std::vector<G4String>  m_typesPhys;

  std::vector<G4String>  m_regnamesSubCut;

  std::vector<G4String>  m_procBiasedXS;
  std::vector<G4double>  m_factBiasedXS;
  std::vector<G4bool>    m_weightBiasedXS;

  std::vector<G4String>  m_procForced;
  std::vector<G4String>  m_regnamesForced;
  std::vector<G4double>  m_lengthForced;
  std::vector<G4bool>    m_weightForced;

  std::vector<G4String>  m_procBiasedSec;
  std::vector<G4String>  m_regnamesBiasedSec;
  std::vector<G4double>  m_factBiasedSec;
  std::vector<G4double>  m_elimBiasedSec;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
