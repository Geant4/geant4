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

G4GEMChannelVI::G4GEMChannelVI(G4int theA, G4int theZ)
  : A(theA), Z(theZ)
{ 
  G4NuclearLevelData* nData = G4NuclearLevelData::GetInstance();
  pairingCorrection = nData->GetPairingCorrection();
  const G4LevelManager* lManager = nullptr;
  if(A > 4) { lManager = nData->GetLevelManager(Z, A); }
  evapMass = G4NucleiProperties::GetNuclearMass(A, Z);
  evapMass2 = evapMass*evapMass;

  cBarrier = new G4CoulombBarrier(A, Z);
  fProbability = new G4GEMProbabilityVI(A, Z, lManager);

  resA = resZ = fragZ = fragA = 0;
  mass = resMass = 0.0; 
}

G4GEMChannelVI::~G4GEMChannelVI()
{
  delete cBarrier;
  delete fProbability;
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

  const G4double exc = fragment->GetExcitationEnergy();
  const G4double delta0 = 
    std::max(pairingCorrection->GetPairingCorrection(fragA, fragZ),0.0);
  if(exc < delta0) { return 0.0; }
 
  resMass = G4NucleiProperties::GetNuclearMass(resA, resZ);
  const G4double fragM = fragment->GetGroundStateMass() + exc;
  const G4double CB = cBarrier->GetCoulombBarrier(resA, resZ, exc);
  
  const G4double delta1 = 
    std::max(0.0,pairingCorrection->GetPairingCorrection(resA,resZ));
  if(fragM <= resMass + CB + delta1) { return 0.0; }

  fProbability->SetDecayKinematics(resZ, resA, resMass, fragM);
  G4double prob = fProbability->ComputeTotalProbability(*fragment, CB);
  //G4cout<<"G4EvaporationChannel: probability= "<< prob <<G4endl;
  return prob;
}

G4Fragment* G4GEMChannelVI::EmittedFragment(G4Fragment* theNucleus)
{
  // assumed, that TotalProbability(...) was already called
  // if value iz zero no possiblity to sample final state
  G4Fragment* evFragment = nullptr;
  G4LorentzVector lv0 = theNucleus->GetMomentum();
  if(resA <= 4 || fProbability->GetProbability() == 0.0) {
    G4double ekin = 
      std::max(0.5*(mass*mass - resMass*resMass + evapMass2)/mass 
               - evapMass, 0.0);
    G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*evapMass))
                       *G4RandomDirection(), ekin + evapMass);
    lv.boost(lv0.boostVector());
    evFragment = new G4Fragment(A, Z, lv);
    lv0 -= lv;
  } else {
    evFragment = fProbability->SampleEvaporationFragment();
    G4LorentzVector lv = evFragment->GetMomentum();
    lv.boost(lv0.boostVector());
    evFragment->SetMomentum(lv);
    lv0 -= lv; 
  }
  theNucleus->SetZandA_asInt(resZ, resA);
  theNucleus->SetMomentum(lv0);

  return evFragment;  
} 

void G4GEMChannelVI::Dump() const
{}



