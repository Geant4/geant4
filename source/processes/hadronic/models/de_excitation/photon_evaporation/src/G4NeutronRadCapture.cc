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
// $Id: G4NeutronRadCapture.cc,v 1.6 2010-11-17 16:21:32 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
#include "G4ParticleDefinition.hh"
#include "G4Fragment.hh"
#include "G4FragmentVector.hh"
#include "G4NucleiProperties.hh"
#include "G4PhotonEvaporation.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"

G4NeutronRadCapture::G4NeutronRadCapture() 
  : G4HadronicInteraction("nRadCapture")
{
  lowestEnergyLimit = 0.1*eV;
  SetMinEnergy( 0.0*GeV );
  SetMaxEnergy( 100.*TeV );
  photonEvaporation = new G4PhotonEvaporation();
  //photonEvaporation = 0;
}

G4NeutronRadCapture::~G4NeutronRadCapture()
{
  delete photonEvaporation;
}

G4HadFinalState* G4NeutronRadCapture::ApplyYourself(
		 const G4HadProjectile& aTrack, G4Nucleus& theNucleus)
{
  theParticleChange.Clear();

  G4int A = theNucleus.GetA_asInt();
  G4int Z = theNucleus.GetZ_asInt();

  // Create initial state
  G4double m1 = G4NucleiProperties::GetNuclearMass(A, Z);
  G4LorentzVector lv0(0.0,0.0,0.0,m1);   
  G4LorentzVector lv1 = aTrack.Get4Momentum() + lv0;

  // simplified method of 1 gamma emission
  if(A <= 3) {

    G4ThreeVector bst = lv1.boostVector();
    G4double M  = lv1.mag();

    ++A;
    G4double m2 = G4NucleiProperties::GetNuclearMass(A, Z);
    if(M - m2 <= lowestEnergyLimit) {
      return &theParticleChange;
    }
 
    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCapture::DoIt: Eini(MeV)=" 
	     << aTrack.GetKineticEnergy()/MeV << "  Eexc(MeV)= " 
	     << (M - m2)/MeV 
	     << "  Z= " << Z << "  A= " << A << G4endl;
    }
    G4double e1 = (M - m2)*(M + m2)/(2*M);
    G4double cost = 2.0*G4UniformRand() - 1.0;
    if(cost > 1.0) {cost = 1.0;}
    else if(cost < -1.0) {cost = -1.0;}
    G4double sint = std::sqrt((1. - cost)*(1.0 + cost));
    G4double phi  = G4UniformRand()*CLHEP::twopi;
    G4LorentzVector lv2(e1*sint*std::cos(phi),e1*sint*std::sin(phi),e1*cost,e1);
    lv2.boost(bst);
    theParticleChange.AddSecondary(new G4DynamicParticle(G4Gamma::Gamma(), lv2));
    G4ParticleDefinition* theDef = 0;

    lv1 -= lv2; 
    if      (Z == 1 && A == 2) {theDef = G4Deuteron::Deuteron();}
    else if (Z == 1 && A == 3) {theDef = G4Triton::Triton();}
    else if (Z == 2 && A == 3) {theDef = G4He3::He3();}
    else if (Z == 2 && A == 4) {theDef = G4Alpha::Alpha();}
    else  
      {
	theDef = 
	  G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
      }

    if (verboseLevel > 1) {
      G4cout << "Gamma 4-mom: " << lv2 << "   " 
	     << theDef->GetParticleName() << "   " << lv1 << G4endl;
    }
    if(theDef) {
      theParticleChange.AddSecondary(new G4DynamicParticle(theDef, lv1));
    }
 
  // Use photon evaporation  
  } else {
 
    G4Fragment aFragment(A+1, Z, lv1);

    if (verboseLevel > 1) {
      G4cout << "G4NeutronRadCapture::ApplyYourself initial G4Fragmet:" << G4endl;
      G4cout << aFragment << G4endl;
    }

    //
    // Sample final state
    //
    G4FragmentVector* fv = photonEvaporation->BreakUpFragment(&aFragment);
    if(fv) {
      size_t n = fv->size();

      if (verboseLevel > 1) {
	G4cout << "G4NeutronRadCapture: " << n << " final particle" << G4endl;
      }
      for(size_t i=0; i<n; ++i) {
	G4Fragment* f = (*fv)[i];    
	G4LorentzVector lvres = f->GetMomentum();   
	Z = f->GetZ_asInt();
	A = f->GetA_asInt();

	G4ParticleDefinition* theDef = 0;
	if(0 == Z && 0 == A) {theDef =  f->GetParticleDefinition();}
	else
	  {
	    theDef = 
	      G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Z,A,0.0);
	  }

	if (verboseLevel > 1) {
	  G4cout << i << ". " << theDef->GetParticleName()
		 << "   " << lvres << G4endl;
	}
	if(theDef) {
	  theParticleChange.AddSecondary(new G4DynamicParticle(theDef, lvres));
	}
	delete f;
      }
      delete fv;
    }
  }
  return &theParticleChange;
}

