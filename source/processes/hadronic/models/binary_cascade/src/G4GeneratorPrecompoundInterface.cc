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
// $Id: G4GeneratorPrecompoundInterface.cc,v 1.11 2010-11-10 17:04:35 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -----------------------------------------------------------------------------
//      GEANT 4 class file
//
//      History: first implementation
//      HPW, 10DEC 98, the decay part originally written by Gunter Folger 
//                in his FTF-test-program.
//
//      M.Kelsey, 28 Jul 2011 -- Replace loop to decay input secondaries
//		with new utility class, simplify cleanup loops
// -----------------------------------------------------------------------------

#include "G4GeneratorPrecompoundInterface.hh"
#include "G4DynamicParticleVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4V3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4FragmentVector.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4DecayKineticTracks.hh"
#include <algorithm>
#include <vector>


G4GeneratorPrecompoundInterface::G4GeneratorPrecompoundInterface(G4VPreCompoundModel* p) 
  : CaptureThreshold(10*MeV)
{
  proton = G4Proton::Proton();
  neutron = G4Neutron::Neutron();
  if(p) { SetDeExcitation(p); }
  else  { SetDeExcitation(new G4PreCompoundModel(new G4ExcitationHandler())); }
}
         
G4GeneratorPrecompoundInterface::~G4GeneratorPrecompoundInterface()
{}


  // choose to calculate excitation energy from energy balance
#define exactExcitationEnergy

G4ReactionProductVector* G4GeneratorPrecompoundInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
	G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;

	// decay the strong resonances
	G4DecayKineticTracks decay(theSecondaries);

	// prepare the fragment
	G4int anA=theNucleus->GetMassNumber();
	G4int aZ=theNucleus->GetCharge();
	G4int numberOfEx = 0;
	G4int numberOfCh = 0;
	G4int numberOfHoles = 0;
	G4double exEnergy = 0.0;
	G4double R = theNucleus->GetNuclearRadius();
	G4ThreeVector exciton3Momentum(0.,0.,0.);

	// loop over secondaries
	unsigned int amax = theSecondaries->size();
#ifdef exactExcitationEnergy
	G4LorentzVector secondary4Momemtum(0,0,0,0);
#endif  
	for(unsigned int list=0; list<amax; ++list)
	{
		G4KineticTrack *aTrack = (*theSecondaries)[list];
		G4ParticleDefinition* part = aTrack->GetDefinition();
		G4double e = aTrack->Get4Momentum().e();
		G4double mass = aTrack->Get4Momentum().mag();
		G4ThreeVector mom = aTrack->Get4Momentum().vect();
		if((part != proton && part != neutron) ||
				(e > mass + CaptureThreshold) ||
				(aTrack->GetPosition().mag() > R)) {
			G4ReactionProduct * theNew = new G4ReactionProduct(part);
			theNew->SetMomentum(mom);
			theNew->SetTotalEnergy(e);
			theTotalResult->push_back(theNew);
#ifdef exactExcitationEnergy
			secondary4Momemtum += aTrack->Get4Momentum();
#endif  
		} else {
			// within the nucleus, neutron or proton
			// now calculate  A, Z of the fragment, momentum, number of exciton states
			++anA;
			++numberOfEx;
			G4int Z = G4int(part->GetPDGCharge()/eplus + 0.1);
			aZ += Z;
			numberOfCh += Z;
			exciton3Momentum += mom;
			exEnergy += (e - mass);
		}
		delete aTrack;
	}
	delete theSecondaries;

	// loop over wounded nucleus
	G4Nucleon * theCurrentNucleon =
			theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
	while(0 != theCurrentNucleon) {
		if(theCurrentNucleon->AreYouHit()) {
			++numberOfHoles;
			++numberOfEx;
			--anA;
			aZ -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/eplus + 0.1);
			exciton3Momentum -= theCurrentNucleon->Get4Momentum().vect();
			exEnergy += theCurrentNucleon->GetBindingEnergy();
		}
		theCurrentNucleon = theNucleus->GetNextNucleon();
	}

	if(0!=anA && 0!=aZ) {
		G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);
#ifdef exactExcitationEnergy        
		// recalculate exEnergy from Energy balance....
		const G4HadProjectile * primary = GetPrimaryProjectile();
		G4double Einitial= primary->Get4Momentum().e()
            		+ G4NucleiProperties::GetNuclearMass(theNucleus->GetMassNumber(),theNucleus->GetCharge());
		G4double Efinal = fMass + secondary4Momemtum.e();
		if ( (Einitial - Efinal) > 0 ) {
				// G4cout << "G4GPI::Propagate() : positive exact excitation Energy "
				//  << (Einitial - Efinal)/MeV << " MeV, exciton estimate " << exEnergy/MeV << " MeV" << G4endl;
			exEnergy=Einitial - Efinal;
		}
		else {
				//  G4cout << "G4GeneratorPrecompoundInterface::Propagate() : negative exact excitation Energy "
				//  << (Einitial - Efinal)/MeV << " MeV, using exciton estimate " << exEnergy/MeV << " MeV" << G4endl;
			exEnergy=0.;
		}
#endif 
		fMass += exEnergy;

#ifdef exactExcitationEnergy        
		G4LorentzVector exciton4Momentum(exciton3Momentum, fMass);
#else
		G4LorentzVector exciton4Momentum(exciton3Momentum,
				std::sqrt(exciton3Momentum.mag2() + fMass*fMass));
#endif    
		if ( exEnergy > 0.0 ) {  // Need to de-excite the remnant nucleus only if excitation energy > 0.
			G4Fragment anInitialState(anA, aZ, exciton4Momentum);
			anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
			anInitialState.SetNumberOfCharged(numberOfCh);
			anInitialState.SetNumberOfHoles(numberOfHoles);

			G4ReactionProductVector * aPreResult = theDeExcitation->DeExcite(anInitialState);

			// fill pre-compound part into the result, and return
			unsigned int amax = aPreResult->size();
			for(unsigned int ll=0; ll<amax; ++ll) {
				theTotalResult->push_back(aPreResult->operator[](ll));
			}
			delete aPreResult;
		} else {  // No excitation energy, we only need to create the remnant nucleus
			G4ParticleDefinition* theKindOfFragment = 0;
			if (anA == 1 && aZ == 0) {
				theKindOfFragment = G4Neutron::NeutronDefinition();
			} else if (anA == 1 && aZ == 1) {
				theKindOfFragment = G4Proton::ProtonDefinition();
			} else if (anA == 2 && aZ == 1) {
				theKindOfFragment = G4Deuteron::DeuteronDefinition();
			} else if (anA == 3 && aZ == 1) {
				theKindOfFragment = G4Triton::TritonDefinition();
			} else if (anA == 3 && aZ == 2) {
				theKindOfFragment = G4He3::He3Definition();
			} else if (anA == 4 && aZ == 2) {
				theKindOfFragment = G4Alpha::AlphaDefinition();;
			} else {
				theKindOfFragment =
						G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(aZ,anA,0.0);
			}
			if (theKindOfFragment != 0) {
				G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
				theNew->SetMomentum(exciton3Momentum);
				theNew->SetTotalEnergy(fMass);
				//theNew->SetFormationTime(??0.??);
				theTotalResult->push_back(theNew);
			}
		}
	}

	return theTotalResult;
}

G4HadFinalState* G4GeneratorPrecompoundInterface::
ApplyYourself(const G4HadProjectile &, G4Nucleus & )
{
  G4cout << "G4GeneratorPrecompoundInterface: ApplyYourself interface called stand-allone."
	 << G4endl;
  G4cout << "This class is only a mediator between generator and precompound"<<G4endl;
  G4cout << "Please remove from your physics list."<<G4endl;
  throw G4HadronicException(__FILE__, __LINE__, "SEVERE: G4GeneratorPrecompoundInterface model interface called stand-allone.");
  return new G4HadFinalState;
}
