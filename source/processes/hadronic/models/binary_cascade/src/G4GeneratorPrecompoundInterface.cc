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
//
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
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"

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
         
G4ReactionProductVector* G4GeneratorPrecompoundInterface::
Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
{
  G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;

  // decay the strong resonances
  G4KineticTrackVector *result1, *secondaries, *result;
  result1=theSecondaries;
  result=new G4KineticTrackVector();
  //G4cout << "### G4GeneratorPrecompoundInterface::Propagate " 
  //	 << result1->size() << " tracks " << theDeExcitation << G4endl;
  for (unsigned int aResult=0; aResult < result1->size(); ++aResult)
    {
      G4ParticleDefinition * pdef;
      pdef=result1->operator[](aResult)->GetDefinition();
      secondaries=0;
      if ( pdef->IsShortLived() )
	{
	  secondaries = result1->operator[](aResult)->Decay();
	}
      if ( 0 == secondaries )
	{
	  result->push_back(result1->operator[](aResult));
	  result1->operator[](aResult)=NULL;	//protect for clearAndDestroy 
	} 
      else
	{
	  unsigned int amax = secondaries->size();
	  for (unsigned int aSecondary=0; aSecondary<amax; ++aSecondary)
	    {
	      result1->push_back(secondaries->operator[](aSecondary));
	    }
	  delete secondaries;
	}
    }
  //G4cout << "Delete tracks" << G4endl;
  std::for_each(result1->begin(), result1->end(), DeleteKineticTrack());
  delete result1;
     
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
  unsigned int amax = result->size();
  for(unsigned int list=0; list<amax; ++list)
    {
      G4KineticTrack *aTrack = result->operator[](list);
      G4ParticleDefinition* part = aTrack->GetDefinition();
      G4double e = aTrack->Get4Momentum().e();
      G4double mass = aTrack->Get4Momentum().mag();
      G4ThreeVector mom = aTrack->Get4Momentum().vect();
      if((part != proton && part != neutron) ||
         (e > mass + CaptureThreshold) ||
	 (aTrack->GetPosition().mag() > R))
	{
	  G4ReactionProduct * theNew = new G4ReactionProduct(part);
	  theNew->SetMomentum(mom);
	  theNew->SetTotalEnergy(e);
	  theTotalResult->push_back(theNew);            
	}
      else
	{
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
    }
     
  // loop over wounded nucleus
  G4Nucleon * theCurrentNucleon = 
    theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : 0;
  while(0 != theCurrentNucleon)
    {
      if(theCurrentNucleon->AreYouHit()) 
	{
	  ++numberOfHoles;
	  ++numberOfEx;
	  --anA;
	  aZ -= G4int(theCurrentNucleon->GetDefinition()->GetPDGCharge()/eplus + 0.1);
	  exciton3Momentum -= theCurrentNucleon->Get4Momentum().vect();
	  exEnergy += theCurrentNucleon->GetBindingEnergy();
	}
      theCurrentNucleon = theNucleus->GetNextNucleon();
    }   
 
  if(0!=anA && 0!=aZ)
    {
      G4double fMass =  G4NucleiProperties::GetNuclearMass(anA, aZ);
      fMass += exEnergy;

      G4LorentzVector exciton4Momentum(exciton3Momentum, 
				       std::sqrt(exciton3Momentum.mag2() + fMass*fMass));
    
      G4Fragment anInitialState(anA, aZ, exciton4Momentum);
      anInitialState.SetNumberOfParticles(numberOfEx-numberOfHoles);
      anInitialState.SetNumberOfCharged(numberOfCh);
      anInitialState.SetNumberOfHoles(numberOfHoles);
      G4ReactionProductVector * aPreResult = theDeExcitation->DeExcite(anInitialState);

      // fill pre-compound part into the result, and return
      unsigned int amax = aPreResult->size();
       for(unsigned int ll=0; ll<amax; ++ll)
	 {
	   theTotalResult->push_back(aPreResult->operator[](ll));
	 }
       delete aPreResult;
    }
     
  std::for_each(result->begin(), result->end(), DeleteKineticTrack());
  delete result;
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
