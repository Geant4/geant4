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
//J.M. Quesada (August2008). Based on:
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//
// Modified:
// 03-09-2008 J.M. Quesada for external choice of inverse cross section option
// 06-09-2008 J.M. Quesada Also external choices have been added for superimposed 
//                 Coulomb barrier (if useSICB is set true, by default is false) 
// 17-11-2010 V.Ivanchenko in constructor replace G4VEmissionProbability by 
//            G4EvaporationProbability and do not new and delete probability
//            object at each call; use G4Pow

#include "G4EvaporationChannel.hh"
#include "G4EvaporationProbability.hh"
#include "G4CoulombBarrier.hh"
#include "G4NuclearLevelData.hh"
#include "G4NucleiProperties.hh"
#include "G4Pow.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4PhysicsModelCatalog.hh"

G4EvaporationChannel::G4EvaporationChannel(G4int anA, G4int aZ, 
					   G4EvaporationProbability* aprob):
  G4VEvaporationChannel(),
  theProbability(aprob),
  theCoulombBarrier(new G4CoulombBarrier(anA, aZ)),
  theA(anA), theZ(aZ)
{ 
  secID = G4PhysicsModelCatalog::GetModelID("model_G4EvaporationChannel");
  evapMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
  evapMass2 = evapMass*evapMass;
  theLevelData = G4NuclearLevelData::GetInstance();
}

G4EvaporationChannel::~G4EvaporationChannel()
{
  delete theCoulombBarrier;
}

void G4EvaporationChannel::Initialise()
{
  theProbability->Initialise();
  G4VEvaporationChannel::Initialise();  
}

G4double G4EvaporationChannel::GetEmissionProbability(G4Fragment* fragment)
{
  theProbability->ResetProbability();
  G4int fragA = fragment->GetA_asInt();
  G4int fragZ = fragment->GetZ_asInt();
  resA = fragA - theA;
  resZ = fragZ - theZ;

  // Only channels which are physically allowed are taken into account 
  if(resA < theA || resA < resZ || resZ < 0 || (resA == theA && resZ < theZ)
     || ((resA > 1) && (resA == resZ || resZ == 0)))
    { return 0.0; }

  G4double exEnergy = fragment->GetExcitationEnergy();
  G4double delta0 = theLevelData->GetPairingCorrection(fragZ,fragA);
  /*  
  G4cout << "G4EvaporationChannel::Initialize Z= "<<theZ<<" A= "<<theA 
  	 << " FragZ= " << fragZ << " FragA= " << fragA 
	 << " exEnergy= " << exEnergy << " d0= " << delta0 << G4endl;
  */
  if(exEnergy < delta0) { return 0.0; }

  G4double fragMass = fragment->GetGroundStateMass();
  mass = fragMass + exEnergy;

  resMass = G4NucleiProperties::GetNuclearMass(resA, resZ);
  ekinmax = 0.5*((mass-resMass)*(mass+resMass) + evapMass2)/mass - evapMass;

  G4double elim = 0.0;
  if(theZ > 0) {
    bCoulomb = theCoulombBarrier->GetCoulombBarrier(resA, resZ, 0.0);

    // for OPTxs >0 penetration under the barrier is taken into account
    elim = (0 != OPTxs) ? bCoulomb*0.6 : bCoulomb;
  }
  /*
  G4cout << "exEnergy= " << exEnergy << " Ec= " << bCoulomb
         << " elim=" << elim << " d0= " << delta0  
	 << " Free= " << mass - resMass - evapMass 
	 << G4endl;
  */
  if(mass <= resMass + evapMass + elim) { return 0.0; }

  G4double ekinmin = 0.0;
  if(elim > 0.0) {
    G4double resM = mass - evapMass - elim;
    ekinmin = 
      std::max(0.5*((mass-resM)*(mass+resM) + evapMass2)/mass - evapMass, 0.0);
  }
  /*  
  G4cout << "Emin= " <<ekinmin<<" Emax= "<<ekinmax
	 << " mass= " << mass << " resM= " << resMass 
	 << " evapM= " << evapMass << G4endl;
  */
  if(ekinmax <= ekinmin) { return 0.0; }

  theProbability->SetDecayKinematics(resZ, resA, resMass, mass);
  G4double prob = theProbability->TotalProbability(*fragment, ekinmin,
                                                   ekinmax, bCoulomb,
                                                   exEnergy - delta0);
  return prob;
}

G4Fragment* G4EvaporationChannel::EmittedFragment(G4Fragment* theNucleus)
{
  G4double ekin = ekinmax;
  // assumed, that TotalProbability(...) was already called
  // if value iz zero no possiblity to sample final state
  if(resA > 4 && theProbability->GetProbability() > 0.0) {
    ekin = theProbability->SampleEnergy();
  }
  ekin = std::max(ekin, 0.0);
  G4LorentzVector lv0 = theNucleus->GetMomentum();
  G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*evapMass))*G4RandomDirection(), 
                     ekin + evapMass);
  lv.boost(lv0.boostVector());

  G4Fragment* evFragment = new G4Fragment(theA, theZ, lv);
  evFragment->SetCreatorModelID(secID);
  lv0 -= lv;
  theNucleus->SetZAandMomentum(lv0, resZ, resA);
  theNucleus->SetCreatorModelID(secID);
  return evFragment; 
} 

G4double G4EvaporationChannel::ComputeInverseXSection(G4Fragment* frag,
                                                      G4double kinEnergy)
{
  ComputeProbability(frag, kinEnergy);
  return theProbability->CrossSection(kinEnergy, bCoulomb);
}

G4double G4EvaporationChannel::ComputeProbability(G4Fragment* frag,
                                                  G4double kinEnergy)
{
  G4int fragA = frag->GetA_asInt();
  G4int fragZ = frag->GetZ_asInt();
  resA = fragA - theA;
  resZ = fragZ - theZ;
  bCoulomb = theCoulombBarrier->GetCoulombBarrier(resA, resZ, 0.0);
  return theProbability->ComputeProbability(kinEnergy, bCoulomb);
}
