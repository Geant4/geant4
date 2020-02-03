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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
//
// J.M.Quesada (14 february 2009) bug fixed in fission width: missing 
//                                parenthesis in the denominator 


#include "G4FissionProbability.hh"
#include "G4PhysicalConstants.hh"
#include "G4NuclearLevelData.hh"
#include "G4PairingCorrection.hh"
#include "G4EvaporationLevelDensityParameter.hh"
#include "G4FissionLevelDensityParameter.hh"
#include "G4Exp.hh"

G4FissionProbability::G4FissionProbability() :
  G4VEmissionProbability(0, 0),
  theEvapLDP(new G4EvaporationLevelDensityParameter),
  theFissLDP(new G4FissionLevelDensityParameter),
  ownEvapLDP(true),
  ownFissLDP(true)
{
  fPairCorr = G4NuclearLevelData::GetInstance()->GetPairingCorrection();
}

G4FissionProbability::~G4FissionProbability()
{
  if (ownEvapLDP) delete theEvapLDP;
  if (ownFissLDP) delete theFissLDP;
}

G4double 
G4FissionProbability::EmissionProbability(const G4Fragment & fragment, 
					  G4double MaximalKineticEnergy)
  // Compute integrated probability of fission channel
{
  if (MaximalKineticEnergy <= 0.0) { return 0.0; }
  G4int A = fragment.GetA_asInt();
  G4int Z = fragment.GetZ_asInt();
  G4double U = fragment.GetExcitationEnergy();
  
  G4double Ucompound = U - fPairCorr->GetPairingCorrection(A,Z);
  G4double Ufission = U - fPairCorr->GetFissionPairingCorrection(A,Z);
  if(Ucompound < 0.0 || Ufission < 0.0) { return 0.0; }
  
  G4double SystemEntropy = 
    2.0*std::sqrt(theEvapLDP->LevelDensityParameter(A,Z,Ucompound)*Ucompound);
	
  G4double afission = theFissLDP->LevelDensityParameter(A,Z,Ufission);

  G4double Cf = 2.0*std::sqrt(afission*MaximalKineticEnergy);
  G4double Exp1 = (SystemEntropy <= 160.0) ? G4Exp(-SystemEntropy) : 0.0; 
  G4double Exp2 = (SystemEntropy-Cf <= 160.0) ? G4Exp(-SystemEntropy+Cf) : 0.0; 

  // JMQ 14/02/09 BUG fixed in fission probability (missing parenthesis 
  // at denominator)
  G4double probability = (Exp1 + (Cf-1.0)*Exp2) / (4.0*pi*afission);

  return probability;
}

