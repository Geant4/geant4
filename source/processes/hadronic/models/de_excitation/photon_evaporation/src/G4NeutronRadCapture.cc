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
// $Id: G4NeutronRadCapture.cc 88191 2015-02-02 17:27:37Z gcosmo $
//
//
// Physics model class G4NeutronRadCapture 
// Created:  31 August 2009
// Author  V.Ivanchenko
//  
// Modified:
// 09.09.2010 V.Ivanchenko added usage of G4PhotonEvaporation 
//

#include "G4NeutronRadCapture.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4VEvaporationChannel.hh"
#include "G4PhotonEvaporation.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Electron.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"

G4NeutronRadCapture::G4NeutronRadCapture() 
  : G4HadronicInteraction("nRadCapture"),
    photonEvaporation(nullptr),lab4mom(0.,0.,0.,0.)
{
  lowestEnergyLimit = 10*CLHEP::eV;
  minExcitation = 0.1*CLHEP::keV;
  SetMinEnergy( 0.0*CLHEP::GeV );
  SetMaxEnergy( 100.*CLHEP::TeV );

  electron = G4Electron::Electron();
  icID = -1;

  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
}

G4NeutronRadCapture::~G4NeutronRadCapture()
{
  delete photonEvaporation;
}

void G4NeutronRadCapture::InitialiseModel()
{
  if(photonEvaporation != nullptr) { return; }
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  minExcitation = param->GetMinExcitation();
  icID = param->GetInternalConversionID();

  photonEvaporation = new G4PhotonEvaporation();
  photonEvaporation->Initialise();
  photonEvaporation->SetICM(true);
}

G4HadFinalState* G4NeutronRadCapture::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  theParticleChange.Clear();
  theParticleChange.SetStatusChange(stopAndKill);

  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();

  G4double time = aTrack.GetGlobalTime();

  // Create initial state
  lab4mom.set(0.,0.,0.,G4NucleiProperties::GetNuclearMass(A, Z));
  lab4mom += aTrack.Get4Momentum();

  G4double M  = lab4mom.mag();
  ++A;
  G4double mass = G4NucleiProperties::GetNuclearMass(A, Z);
  //G4cout << "Capture start: Z= " << Z << " A= " << A 
  //	 << " LabM= " << M << " Mcompound= " << mass << G4endl;

  // simplified method of 1 gamma emission
  if(A <= 4) {

    G4ThreeVector bst = lab4mom.boostVector();

    if(M - mass <= lowestEnergyLimit) {
      return &theParticleChange;
    }
 
    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCapture::DoIt: Eini(MeV)=" 
	     << aTrack.GetKineticEnergy()/MeV << "  Eexc(MeV)= " 
	     << (M - mass)/MeV 
	     << "  Z= " << Z << "  A= " << A << G4endl;
    }
    G4double e1 = (M - mass)*(M + mass)/(2*M);

    G4double cost = 2.0*G4UniformRand() - 1.0;
    if(cost > 1.0) {cost = 1.0;}
    else if(cost < -1.0) {cost = -1.0;}
    G4double sint = std::sqrt((1. - cost)*(1.0 + cost));
    G4double phi  = G4UniformRand()*CLHEP::twopi;

    G4LorentzVector lv2(e1*sint*std::cos(phi),e1*sint*std::sin(phi),
			e1*cost,e1);
    lv2.boost(bst);
    G4HadSecondary* news = 
      new G4HadSecondary(new G4DynamicParticle(G4Gamma::Gamma(), lv2));
    news->SetTime(time);
    theParticleChange.AddSecondary(*news);
    delete news;

    const G4ParticleDefinition* theDef = 0;

    lab4mom -= lv2; 
    if      (Z == 1 && A == 2) {theDef = G4Deuteron::Deuteron();}
    else if (Z == 1 && A == 3) {theDef = G4Triton::Triton();}
    else if (Z == 2 && A == 3) {theDef = G4He3::He3();}
    else if (Z == 2 && A == 4) {theDef = G4Alpha::Alpha();}
    else {  theDef = theTableOfIons->GetIon(Z,A,0.0,noFloat,0); }

    if (verboseLevel > 1) {
      G4cout << "Gamma 4-mom: " << lv2 << "   " 
	     << theDef->GetParticleName() << "   " << lab4mom << G4endl;
    }
    if(theDef) {
      news = new G4HadSecondary(new G4DynamicParticle(theDef, lab4mom));
      news->SetTime(time);
      theParticleChange.AddSecondary(*news);
      delete news;
    }
 
  // Use photon evaporation  
  } else {
 
    // protection against wrong kinematic 
    if(M < mass) {
      G4double etot = std::max(mass, lab4mom.e());
      G4double ptot = std::sqrt((etot - mass)*(etot + mass));
      G4ThreeVector v = lab4mom.vect().unit();
      lab4mom.set(v.x()*ptot,v.y()*ptot,v.z()*ptot,etot);
    }

    G4Fragment* aFragment = new G4Fragment(A, Z, lab4mom);

    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCapture::ApplyYourself initial G4Fragmet:" 
	     << G4endl;
      G4cout << aFragment << G4endl;
    }

    //
    // Sample final state
    //
    G4FragmentVector* fv = photonEvaporation->BreakUpFragment(aFragment);
    if(!fv) { fv = new G4FragmentVector(); }
    fv->push_back(aFragment);
    size_t n = fv->size();

    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCapture: " << n << " final particle icID= " << icID << G4endl;
    }
    for(size_t i=0; i<n; ++i) {

      G4Fragment* f = (*fv)[i];    
      G4double etot = f->GetMomentum().e();

      Z = f->GetZ_asInt();
      A = f->GetA_asInt();

      const G4ParticleDefinition* theDef;
      if(0 == Z && 0 == A) {theDef =  f->GetParticleDefinition();}
      else if (Z == 1 && A == 2) {theDef = G4Deuteron::Deuteron();}
      else if (Z == 1 && A == 3) {theDef = G4Triton::Triton();}
      else if (Z == 2 && A == 3) {theDef = G4He3::He3();}
      else if (Z == 2 && A == 4) {theDef = G4Alpha::Alpha();}
      else {
        G4double eexc = f->GetExcitationEnergy();
	if(eexc <= minExcitation) { eexc = 0.0; }
	theDef = theTableOfIons->GetIon(Z, A, eexc, noFloat, 0);
	/*	
	G4cout << "### NC Find ion Z= " << Z << " A= " << A
	       << " Eexc(MeV)= " << eexc/MeV << "  " 
	       << theDef << G4endl;
	*/
      }
      G4double ekin = std::max(0.0,etot - theDef->GetPDGMass());
      if (verboseLevel > 1) {
	G4cout << i << ". " << theDef->GetParticleName()
	       << " Ekin(MeV)= " << etot/MeV
	       << " p: " << f->GetMomentum().vect() 
	       << G4endl;
      }
      G4HadSecondary* news = new G4HadSecondary(
	new G4DynamicParticle(theDef,
			      f->GetMomentum().vect().unit(),
			      ekin));
      G4double timeF = f->GetCreationTime();
      if(timeF < 0.0) { timeF = 0.0; }
      news->SetTime(time + timeF);
      if(theDef == electron) { news->SetCreatorModelType(icID); }
      theParticleChange.AddSecondary(*news);
      delete news;
      delete f;
    }
    delete fv;
  }
  //G4cout << "Capture done" << G4endl;
  return &theParticleChange;
}

