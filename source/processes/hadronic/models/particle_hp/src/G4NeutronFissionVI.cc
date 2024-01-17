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
// Physics model class G4NeutronFissionVI
// Created:  03 October 2023
// Author  V.Ivanchenko
//  

#include "G4NeutronFissionVI.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4VPreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4PhotonEvaporation.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleHPManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Electron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"
#include "G4HadronicParameters.hh"
#include "G4PhysicsModelCatalog.hh"

G4NeutronFissionVI::G4NeutronFissionVI() 
  : G4HadronicInteraction("nFissionVI"),
    fManagerHP(G4ParticleHPManager::GetInstance()),
    minExcitation(0.1*CLHEP::keV),
    emaxT(fManagerHP->GetMaxEnergyDoppler()),
    lab4mom(0.,0.,0.,0.)
{
  SetMinEnergy( 0.0*CLHEP::GeV );
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
}

G4NeutronFissionVI::~G4NeutronFissionVI()
{
  if (fLocalHandler) { delete fHandler; }
}

void G4NeutronFissionVI::InitialiseModel()
{
  if (fFission != nullptr && fHandler != nullptr) { return; }
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  if (nullptr != p) {
    fHandler = (static_cast<G4VPreCompoundModel*>(p))->GetExcitationHandler();
  }
  if (nullptr == fHandler) {
    fHandler = new G4ExcitationHandler();
    fLocalHandler = true;
  }
  fHandler->Initialise();
  fFission = fHandler->GetEvaporation()->GetFissionChannel();

  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  minExcitation = param->GetMinExcitation();
  secID = G4PhysicsModelCatalog::GetModelID("model_" + GetModelName());
}

G4HadFinalState* G4NeutronFissionVI::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  theParticleChange.Clear();
  theParticleChange.SetStatusChange(stopAndKill);
  G4double T = aTrack.GetMaterial()->GetTemperature();

  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();

  G4double time = aTrack.GetGlobalTime();
  G4double ekin = aTrack.GetKineticEnergy();

  // Create initial state
  G4double mass = G4NucleiProperties::GetNuclearMass(A, Z);

  // no Doppler broading
  G4double factT = T/CLHEP::STP_Temperature; 

  if (ekin >= emaxT*factT || fManagerHP->GetNeglectDoppler()) {
    lab4mom.set(0.,0.,0.,mass);

  } else {
    G4double lambda = 1.0/(CLHEP::k_Boltzmann*T);
    G4double erand = G4RandGamma::shoot(2.0, lambda);
    auto mom = G4RandomDirection()*std::sqrt(2*mass*erand);
    lab4mom.set(mom.x(), mom.y(), mom.z(), mass + erand);
  }

  lab4mom += aTrack.Get4Momentum();

  G4double M = lab4mom.mag();
  ++A;
  mass = G4NucleiProperties::GetNuclearMass(A, Z);
  //G4cout << "Fission start: Z= " << Z << " A= " << A 
  //	 << " LabM= " << M << " Mcompound= " << mass << G4endl;

  // protection against wrong kinematic 
  if (M < mass) {
    G4double etot = std::max(mass, lab4mom.e());
    G4double ptot = std::sqrt((etot - mass)*(etot + mass));
    G4ThreeVector v = lab4mom.vect().unit();
    lab4mom.set(v.x()*ptot,v.y()*ptot,v.z()*ptot,etot);
  }

  G4Fragment* aFragment = new G4Fragment(A, Z, lab4mom);

  if (verboseLevel > 1) {
    G4cout << "G4NeutronFissionVI::ApplyYourself initial G4Fragmet:" 
	   << G4endl;
    G4cout << aFragment << G4endl;
  }

  //
  // Sample final state
  //
  fFission->GetEmissionProbability(aFragment);
  G4Fragment* frag = fFission->EmittedFragment(aFragment);
  G4ReactionProductVector* final = fHandler->BreakItUp(*aFragment);
  if (nullptr != frag) {
    G4ReactionProductVector* v = fHandler->BreakItUp(*frag);
    for (auto & p : *v) {
      final->push_back(p);
    }
    delete v;
    delete frag;
  }

  if (verboseLevel > 1) {
      G4cout << "G4NeutronFissionVI: " << final->size()
             << " final particle secID= " << secID << G4endl;
  }
  for (auto const & ptr : *final) {
    G4double etot = ptr->GetTotalEnergy();
    const G4ParticleDefinition* theDef = ptr->GetDefinition();
    ekin = std::max(0.0, etot - theDef->GetPDGMass());
    if (verboseLevel > 1) {
      G4cout << theDef->GetParticleName()
	     << " Ekin(MeV)= " << ekin/MeV
	     << " p: " << ptr->GetMomentum() 
	     << G4endl;
    }
    G4HadSecondary* news = new G4HadSecondary(
      new G4DynamicParticle(theDef, ptr->GetMomentum().unit(), ekin));
    G4double timeF = std::max(ptr->GetFormationTime(), 0.0);
    news->SetTime(time + timeF);
    news->SetCreatorModelID(secID);
    theParticleChange.AddSecondary(*news);
    delete news;
    delete ptr;
  }
  delete final;
  delete aFragment;

  //G4cout << "Fission done" << G4endl;
  return &theParticleChange;
}

