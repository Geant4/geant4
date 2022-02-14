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
#include "G4Alpha.hh"
#include "G4PhysicsModelCatalog.hh"

G4EvaporationChannel::G4EvaporationChannel(G4int anA, G4int aZ, 
					   G4EvaporationProbability* aprob):
  G4VEvaporationChannel(),
  theA(anA),
  theZ(aZ),
  secID(-1),
  theProbability(aprob),
  theCoulombBarrier(new G4CoulombBarrier(anA, aZ))
{ 
  resA = resZ = 0;
  secID = G4PhysicsModelCatalog::GetModelID("model_G4EvaporationChannel");
  mass = resMass = 0.0; 
  evapMass = G4NucleiProperties::GetNuclearMass(theA, theZ);
  //G4cout << "G4EvaporationChannel: Z= " << theZ << " A= " << theA 
  //      << " M(GeV)= " << evapMass/GeV << G4endl;
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
  G4double bCoulomb = 0.0;
  G4double elim = 0.0;
  if(theZ > 0) {
    bCoulomb = theCoulombBarrier->GetCoulombBarrier(resA,resZ,exEnergy);

    // for OPTxs >0 penetration under the barrier is taken into account
    const G4double dCB = 3.5*CLHEP::MeV;
    elim = (0 != OPTxs) ? 
      std::max(bCoulomb*0.5, bCoulomb - dCB*theZ) : bCoulomb;
  }
  /*
  G4cout << "exEnergy= " << exEnergy << " Ec= " << bCoulomb
         << " d0= " << delta0 
	 << " Free= " << mass - resMass - evapMass 
	 << G4endl;
  */
  if(mass <= resMass + evapMass + elim) { return 0.0; }

  G4double twoMass = mass + mass;
  G4double ekinmax = 
    ((mass-resMass)*(mass+resMass) + evapMass2)/twoMass - evapMass;
  G4double ekinmin = 0.0;
  if(elim > 0.0) {
    G4double resM = std::max(mass - evapMass - elim, resMass);
    ekinmin = 
      std::max(((mass-resM)*(mass+resM) + evapMass2)/twoMass - evapMass,0.0);
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
  /*  
  G4cout<<"G4EvaporationChannel: prob= "<< prob << " Z= " << theZ 
        << " A= " << theA << " E1= " << ekinmin << " E2= " << ekinmax 
        << G4endl;
  */
  return prob;
}

G4Fragment* G4EvaporationChannel::EmittedFragment(G4Fragment* theNucleus)
{
  G4double ekin;
  // assumed, that TotalProbability(...) was already called
  // if value iz zero no possiblity to sample final state
  if(resA <= 4 || theProbability->GetProbability() == 0.0) {
    ekin = 0.5*(mass*mass - resMass*resMass + evapMass2)/mass - evapMass;
  } else {
    ekin = theProbability->SampleEnergy();
  }
  ekin = std::max(ekin, 0.0);
  G4LorentzVector lv0 = theNucleus->GetMomentum();
  G4LorentzVector lv(std::sqrt(ekin*(ekin + 2.0*evapMass))*G4RandomDirection(), 
                     ekin + evapMass);
  lv.boost(lv0.boostVector());

  G4Fragment* evFragment = new G4Fragment(theA, theZ, lv);
  if(evFragment != nullptr) { evFragment->SetCreatorModelID(secID); }
  lv0 -= lv;
  theNucleus->SetZandA_asInt(resZ, resA);
  theNucleus->SetMomentum(lv0);
  theNucleus->SetCreatorModelID(secID);
  
  //G4cout << "Residual: Z= " << resZ << " A= " << resA << " Eex= " 
  //	 << theNucleus->GetExcitationEnergy() << G4endl;
  return evFragment; 
} 
