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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modifications:
// 28.10.2010 V.Ivanchenko defined members in constructor and cleaned up

#include "G4VEmissionProbability.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4DeexPrecoParameters.hh"
#include "Randomize.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

G4VEmissionProbability::G4VEmissionProbability(G4int Z, G4int A)
  : pVerbose(1), theZ(Z), theA(A)
{
  pNuclearLevelData = G4NuclearLevelData::GetInstance(); 
  pG4pow = G4Pow::GetInstance();
  if(A > 0) { pEvapMass = G4NucleiProperties::GetNuclearMass(theA, theZ); }
  G4DeexPrecoParameters* param = pNuclearLevelData->GetParameters();
  OPTxs = param->GetDeexModelType();
}

void G4VEmissionProbability::Initialise()
{
  G4DeexPrecoParameters* param = pNuclearLevelData->GetParameters();
  pVerbose = param->GetVerbose();
  fFD = param->GetDiscreteExcitationFlag();
  fMaxLifeTime = param->GetMaxLifeTime();
  pTolerance = param->GetMinExcitation();
  pWidth = param->GetNuclearLevelWidth();
}

void G4VEmissionProbability::ResetIntegrator(G4double de, G4double eps)
{
  InitialiseIntegrator(eps, 0.25, 1.25, de, 0.1*CLHEP::MeV, 2*CLHEP::MeV);
}

G4double G4VEmissionProbability::EmissionProbability(const G4Fragment&, G4double)
{
  return 0.0;
}

G4double G4VEmissionProbability::ComputeProbability(G4double, G4double)
{
  return 0.0;
}

G4double G4VEmissionProbability::IntegrateProbability(G4double elow, 
                                                      G4double ehigh, 
                                                      G4double cb)
{
  pProbability = 0.0;
  if (elow >= ehigh) { return pProbability; }

  eCoulomb = cb;
  pProbability = ComputeIntegral(elow, ehigh);

  if (pVerbose > 1) { 
    G4cout << "G4VEmissionProbability::IntegrateProbability Probability="
	   << pProbability << " Z=" << theZ << " A=" << theA << G4endl; 
  }
  return pProbability;
}

G4double G4VEmissionProbability::SampleEnergy()
{
  G4double ekin = SampleValue();
  G4double enew = FindRecoilExcitation(ekin);
  if (pVerbose > 1) {
    G4cout << "### G4VEmissionProbability::SampleEnergy: Efin(MeV)= " 
	   << enew << " E=" << ekin << "  Eexc=" << fExcRes << G4endl;
  }
  return enew;
}

G4double G4VEmissionProbability::ProbabilityDensityFunction(G4double e)
{
  return ComputeProbability(e, eCoulomb);
}

G4double G4VEmissionProbability::FindRecoilExcitation(const G4double e)
{
  G4double mass = pEvapMass + fExc;
    
  G4double m02 = pMass*pMass;
  G4double m12 = mass*mass;
  G4double m22 = pResMass*pResMass;
  G4double mres = std::sqrt(m02 + m12 - 2.*pMass*(mass + e));

  fExcRes = mres - pResMass;

  if (pVerbose > 1) {
    G4cout << "### FindRecoilExcitation for resZ= " 
           << resZ << " resA= " << resA 
           << " evaporated Z= " << theZ << " A= " << theA
	   << " Ekin= " << e << " Eexc= " << fExcRes << G4endl;
  }

  // residual nucleus is in the ground state
  if(fExcRes < pTolerance) {
    fExcRes = 0.0;
    return std::max(0.5*(m02 + m12 - m22)/pMass - mass, 0.0);
  }
  if (!fFD) { return e; }
 
  // select final state excitation
  auto lManager = pNuclearLevelData->GetLevelManager(resZ, resA);
  if(nullptr == lManager) { return e; }

  // levels are not known
  if(fExcRes > lManager->MaxLevelEnergy() + pTolerance) { return e; }

  // find level
  std::size_t idx = lManager->NearestLevelIndex(fExcRes);
  auto level = lManager->GetLevel(idx);
  G4double ltime = level->GetTimeGamma();
  G4double elevel = lManager->LevelEnergy(idx);

  G4double efinal = e;

  // is possible to use level energy?
  if ((idx <= 1 || std::abs(elevel - fExcRes) <= pWidth || ltime >= fMaxLifeTime) &&
      (pMass >= mass + pResMass + elevel)) { 
    G4double massR = pResMass + elevel;
    G4double mr2 = massR*massR;
    fExcRes = elevel;
    efinal = std::max(0.5*(m02 + m12 - mr2)/pMass - mass, 0.0);
  }
  return efinal;
}
