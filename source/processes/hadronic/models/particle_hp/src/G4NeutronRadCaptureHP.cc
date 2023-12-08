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
// Physics model class G4NeutronRadCaptureHP
//         derived from G4NeutronRadCapture
//
// Created:  02 October 2023
// Author  V.Ivanchenko
//  
// 

#include "G4NeutronRadCaptureHP.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4VPreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4VEvaporationChannel.hh"
#include "G4PhotonEvaporation.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleHPManager.hh"
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

G4NeutronRadCaptureHP::G4NeutronRadCaptureHP() 
  : G4HadronicInteraction("nRadCaptureHP"),
    electron(G4Electron::Electron()),
    fManagerHP(G4ParticleHPManager::GetInstance()),
    lowestEnergyLimit(1.0e-11*CLHEP::eV),
    minExcitation(0.1*CLHEP::keV),
    emax(20*CLHEP::MeV),
    emaxT(fManagerHP->GetMaxEnergyDoppler()),
    lab4mom(0.,0.,0.,0.)
{
  SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
}

G4NeutronRadCaptureHP::~G4NeutronRadCaptureHP()
{
  if (fLocalPE) { delete photonEvaporation; }
}

void G4NeutronRadCaptureHP::BuildPhysicsTable(const G4ParticleDefinition&)
{
  if (photonEvaporation != nullptr) { return; }
  G4HadronicInteraction* p =
    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
  if (nullptr != p) {
    auto handler =
      (static_cast<G4VPreCompoundModel*>(p))->GetExcitationHandler();
    if (nullptr != handler)
      photonEvaporation = handler->GetPhotonEvaporation();
  }
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  minExcitation = param->GetMinExcitation();
  icID = G4PhysicsModelCatalog::GetModelID("model_e-InternalConversion");
  secID = G4PhysicsModelCatalog::GetModelID("model_" + GetModelName());
  if (nullptr == photonEvaporation) {
    photonEvaporation = new G4PhotonEvaporation();
    fLocalPE = true;
  }
  photonEvaporation->Initialise();
  photonEvaporation->SetICM(true);
}

G4HadFinalState* G4NeutronRadCaptureHP::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  theParticleChange.Clear();
  G4double ekin = aTrack.GetKineticEnergy();
  if (ekin > emax) {
    return &theParticleChange;
  }

  theParticleChange.SetStatusChange(stopAndKill);
  G4double T = aTrack.GetMaterial()->GetTemperature();

  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();

  G4double time = aTrack.GetGlobalTime();

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
  //G4cout << "Capture start: Z= " << Z << " A= " << A 
  //	 << " LabM= " << M << " Mcompound= " << mass << G4endl;

  // simplified method of 1 gamma emission
  if (A <= 4) {
 
    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCaptureHP::DoIt: Eini(MeV)=" 
	     << ekin/MeV << "  Eexc(MeV)= " 
	     << (M - mass)/MeV 
	     << "  Z= " << Z << "  A= " << A << G4endl;
    }
    if (M - mass > lowestEnergyLimit) {
      G4ThreeVector bst = lab4mom.boostVector();
      G4double e1 = (M - mass)*(M + mass)/(2*M);
      G4LorentzVector lv2(e1*G4RandomDirection(),e1);
      lv2.boost(bst);
      if (verboseLevel > 1) {
        G4cout << "Gamma 4-mom: " << lv2 << " Escm(MeV)=" << e1/CLHEP::MeV << G4endl;
      }
      lab4mom -= lv2; 
      G4HadSecondary* news = 
	new G4HadSecondary(new G4DynamicParticle(G4Gamma::Gamma(), lv2));
      news->SetTime(time);
      news->SetCreatorModelID(secID);
      theParticleChange.AddSecondary(*news);
      delete news;
    }

    const G4ParticleDefinition* theDef = nullptr;

    if      (Z == 1 && A == 2) {theDef = G4Deuteron::Deuteron(); }
    else if (Z == 1 && A == 3) {theDef = G4Triton::Triton(); }
    else if (Z == 2 && A == 3) {theDef = G4He3::He3(); }
    else if (Z == 2 && A == 4) {theDef = G4Alpha::Alpha(); }
    else { theDef = theTableOfIons->GetIon(Z, A, 0.0, noFloat, 0); }

    if (nullptr != theDef) {
      G4HadSecondary* news =
        new G4HadSecondary(new G4DynamicParticle(theDef, lab4mom));
      news->SetTime(time);
      news->SetCreatorModelID(secID);
      theParticleChange.AddSecondary(*news);
      delete news;
    }
 
  // Use photon evaporation  
  } else {
 
    // protection against wrong kinematic 
    if (M < mass) {
      G4double etot = std::max(mass, lab4mom.e());
      G4double ptot = std::sqrt((etot - mass)*(etot + mass));
      G4ThreeVector v = lab4mom.vect().unit();
      lab4mom.set(v.x()*ptot,v.y()*ptot,v.z()*ptot,etot);
    }

    G4Fragment* aFragment = new G4Fragment(A, Z, lab4mom);

    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCaptureHP::ApplyYourself initial G4Fragmet:" 
	     << G4endl;
      G4cout << aFragment << G4endl;
    }

    //
    // Sample final state
    //
    G4FragmentVector* fv = photonEvaporation->BreakUpFragment(aFragment);
    if (nullptr == fv) { fv = new G4FragmentVector(); }
    fv->push_back(aFragment);

    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCaptureHP: " << fv->size() << " final particle icID= "
	     << icID << G4endl;
    }
    for (auto const & f : *fv) {
      G4double etot = f->GetMomentum().e();

      Z = f->GetZ_asInt();
      A = f->GetA_asInt();

      const G4ParticleDefinition* theDef;
      if (0 == Z && 0 == A) { theDef =  f->GetParticleDefinition(); }
      else if (Z == 1 && A == 2) { theDef = G4Deuteron::Deuteron(); }
      else if (Z == 1 && A == 3) { theDef = G4Triton::Triton(); }
      else if (Z == 2 && A == 3) { theDef = G4He3::He3(); }
      else if (Z == 2 && A == 4) { theDef = G4Alpha::Alpha(); }
      else {
        G4double eexc = f->GetExcitationEnergy();
	if (eexc <= minExcitation) { eexc = 0.0; }
	theDef = theTableOfIons->GetIon(Z, A, eexc, noFloat, 0);
	/*	
	G4cout << "### NC Find ion Z= " << Z << " A= " << A
	       << " Eexc(MeV)= " << eexc/MeV << "  " 
	       << theDef << G4endl;
	*/
      }
      ekin = std::max(0.0, etot - theDef->GetPDGMass());
      if (verboseLevel > 1) {
	G4cout << theDef->GetParticleName()
	       << " Ekin(MeV)= " << ekin/MeV
	       << " p: " << f->GetMomentum().vect() 
	       << G4endl;
      }
      G4HadSecondary* news =
        new G4HadSecondary(new G4DynamicParticle(theDef,
                                                 f->GetMomentum().vect().unit(),
                                                 ekin));
      G4double timeF = std::max(f->GetCreationTime(), 0.0);
      news->SetTime(time + timeF);
      if (theDef == electron) { 
        news->SetCreatorModelID(icID); 
      } else {
        news->SetCreatorModelID(secID);
      }
      theParticleChange.AddSecondary(*news);
      delete news;
      delete f;
    }
    delete fv;
  }
  //G4cout << "Capture done" << G4endl;
  return &theParticleChange;
}

