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
// GEM de-excitation model
// by V. Ivanchenko (July 2019)
//

#include "G4GEMChannelVI.hh"
#include "G4GEMProbabilityVI.hh"
#include "G4VCoulombBarrier.hh"
#include "G4CoulombBarrier.hh"
#include "G4PairingCorrection.hh"
#include "G4NuclearLevelData.hh"
#include "G4LevelManager.hh"
#include "G4NucleiProperties.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicsModelCatalog.hh"
#include "Randomize.hh"

namespace
{
  const G4double minExc = 1.0*CLHEP::MeV;
  const G4int nProbMax = 10;
}

G4GEMChannelVI::G4GEMChannelVI(G4int theA, G4int theZ)
  : A(theA), Z(theZ)
{
  G4NuclearLevelData* nData = G4NuclearLevelData::GetInstance();
  pairingCorrection = nData->GetPairingCorrection();
  const G4LevelManager* lManager = nullptr;
  if (A > 4) { lManager = nData->GetLevelManager(Z, A); }
  fEvapMass = G4NucleiProperties::GetNuclearMass(A, Z);
  fEvapMass2 = fEvapMass*fEvapMass;

  cBarrier = new G4CoulombBarrier(A, Z);
  fProbability = new G4GEMProbabilityVI(A, Z, lManager);

  fCoeff = CLHEP::millibarn/((CLHEP::pi*CLHEP::hbarc)*(CLHEP::pi*CLHEP::hbarc)); 

  secID = G4PhysicsModelCatalog::GetModelID("model_G4GEMChannelVI");
  if (Z == 0 && A == 1) {
    indexC = 0;
    fCoeff *= 2.0;
  } else if (Z == 1 && A == 1) {
    indexC = 1;
    fCoeff *= 2.0;
  } else if (Z == 1 && A == 2) {
    indexC = 2;
    fCoeff *= 3.0;
  } else if (Z == 1 && A == 3) {
    indexC = 3;
    fCoeff *= 2.0;
  } else if (Z == 2 && A == 3) {
    indexC = 4;
    fCoeff *= 2.0;
  } else if (Z == 2 && A == 4) {
    indexC = 5;
  } else {
    indexC = 6;
  }
}

G4GEMChannelVI::~G4GEMChannelVI()
{
  delete cBarrier;
  delete fProbability;
}

void G4GEMChannelVI::Initialise()
{
  fProbability->Initialise();
  G4VEvaporationChannel::Initialise(); 
}

G4double G4GEMChannelVI::GetEmissionProbability(G4Fragment* fragment)
{
  fProbability->ResetProbability();
  fragZ = fragment->GetZ_asInt();
  fragA = fragment->GetA_asInt();
  resZ = fragZ - Z;
  resA = fragA - A;
  if(resA < A || resA < resZ || resZ < 0 || (resA == A && resZ < Z)) { 
    return 0.0; 
  }

  fExc = fragment->GetExcitationEnergy();
  fMass = fragment->GetGroundStateMass() + fExc;
  fResMass = G4NucleiProperties::GetNuclearMass(resA, resZ);

  // limit for the case when both evaporation and residual 
  // fragments are in ground states
  if (fMass <= fEvapMass + fResMass) { return 0.0; } 

  if (Z > 0) {
    bCoulomb = cBarrier->GetCoulombBarrier(resA, resZ, 0.0);
  }
  G4double de = fMass - fEvapMass - fResMass - bCoulomb;
  nProb = (G4int)(de/minExc);
  if (nProb <= 1 || indexC < 6 || resA <= 4) {
    nProb = 1;
  } else {
    nProb = std::min(nProb, nProbMax);
  }
  if (2 < fVerbose) {
    G4cout << "## G4GEMChannelVI::GetEmissionProbability fragZ="
	   << fragZ << " fragA=" << fragA << " Z=" << Z << " A=" << A
	   << " Eex(MeV)=" << fExc << " nProb=" << nProb 
	   << G4endl;
  }
  fProbability->SetDecayKinematics(resZ, resA, fResMass, fMass);
  G4double sump = 0.0; 
  for (G4int i=0; i<nProb; ++i) {
    G4double exc = std::min(minExc*i, de);
    G4double m1 = fEvapMass + exc;
    G4double e2 = 0.5*((fMass-fResMass)*(fMass+fResMass) + m1*m1)/fMass - m1;
    G4double m2 = fMass - m1 - 0.5*bCoulomb;
    if (m2 < fResMass) {
      nProb = i;
      break;
    }
    G4double e1 = std::max(0.5*((fMass-m2)*(fMass+m2) + m1*m1)/fMass - m1, 0.0);
    if (e1 >= e2) {
      nProb = i;
      break;
    }
    sump += fProbability->TotalProbability(*fragment, e1, e2, bCoulomb, fExc, exc);
    fEData[i].exc = exc;
    fEData[i].ekin1 = e1;
    fEData[i].ekin2 = e2;
    fEData[i].prob = sump;
  }
  return sump;
}

G4Fragment* G4GEMChannelVI::EmittedFragment(G4Fragment* theNucleus)
{
  // assumed, that TotalProbability(...) was already called
  // if value iz zero no possiblity to sample final state
  G4Fragment* evFragment = nullptr;
  G4LorentzVector lv0 = theNucleus->GetMomentum();
  G4double ekin;
  G4double exc = 0.0;
  G4double probMax = std::max(fEData[nProb - 1].prob, 0.0); 
  if (0.0 >= probMax) {
    ekin = std::max(0.5*(fMass*fMass - fResMass*fResMass + fEvapMass2)
		    /fMass - fEvapMass, 0.0);
  } else if (1 == nProb) {
    ekin = fProbability->SampleEnergy(fEData[0].ekin1, fEData[0].ekin2,
				      bCoulomb, fExc, 0.0);
  } else {
    G4double p = G4UniformRand()*probMax;
    G4int i{1};
    for (; i<nProb; ++i) {
      if (p <= fEData[i].prob) { break; }
    }
    G4double e1 = fEData[i - 1].exc;
    G4double e2 = fEData[i].exc;
    G4double p1 = fEData[i - 1].prob;
    G4double p2 = fEData[i].prob;
    exc = e1 + (e2 - e1)*(p - p1)/(p2 - p1);
    ekin = fProbability->SampleEnergy(fEData[i].ekin1, fEData[i].ekin2,
                                      bCoulomb, fExc, exc);
  }
  G4double m1 = fEvapMass + exc;
  G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*m1))
		     *G4RandomDirection(), ekin + m1);
  lv.boost(lv0.boostVector());
  evFragment = new G4Fragment(A, Z, lv);
  lv0 -= lv;
  evFragment->SetCreatorModelID(secID);
  theNucleus->SetZandA_asInt(resZ, resA);
  theNucleus->SetMomentum(lv0);
  theNucleus->SetCreatorModelID(secID);
  
  return evFragment;  
} 

void G4GEMChannelVI::Dump() const
{}



