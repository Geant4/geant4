#include "G4BinaryLightIonReaction.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include <algorithm>
#include "G4ReactionProductVector.hh"
#include <vector>
#include "G4ping.hh"
#include "G4Delete.hh"
  
  G4BinaryLightIonReaction::G4BinaryLightIonReaction()
  : theModel(), theHandler(), theProjectileFragmentation(&theHandler) {}
  
  G4VParticleChange *G4BinaryLightIonReaction::
  ApplyYourself(const G4Track &aTrack, G4Nucleus & targetNucleus )
  {    
    G4ping debug("debug_G4BinaryLightIonReaction");
    G4double a1=aTrack.GetDefinition()->GetBaryonNumber();
    G4double z1=aTrack.GetDefinition()->GetPDGCharge();
    G4double m1=aTrack.GetDefinition()->GetPDGMass();
    G4double a2=targetNucleus.GetN();
    G4double z2=targetNucleus.GetZ();
    G4double m2=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z2, a2);
    debug.push_back(a1);
    debug.push_back(z1);
    debug.push_back(m1);
    debug.push_back(a2);
    debug.push_back(z2);
    debug.push_back(m2);
    G4LorentzVector mom(aTrack.GetDynamicParticle()->Get4Momentum());
    debug.push_back(mom);
    debug.dump();
    G4LorentzRotation toBreit(mom.boostVector());
        
    G4bool swapped = false;
    if(a2<a1)
    {
      swapped = true;
      G4double tmp(0);
      tmp = a2; a2=a1; a1=tmp;
      tmp = z2; z2=z1; z1=tmp;
      tmp = m2; m2=m1; m1=tmp;
      G4LorentzVector it(m1, G4ThreeVector(0,0,0));
      mom = toBreit*it;
    }
    debug.push_back("After swap");
    debug.push_back(a1);
    debug.push_back(z1);
    debug.push_back(a2);
    debug.push_back(z2);
    debug.push_back(mom);
    debug.dump();

    G4ReactionProductVector * result = 0;
    G4V3DNucleus * fancyNucleus(0);
    while(!result)
    {
      G4V3DNucleus * projectile = new G4Fancy3DNucleus;
      projectile->Init(a1, z1);
      fancyNucleus = new G4Fancy3DNucleus;  
      fancyNucleus->Init(a2, z2);

      G4double impactMax = fancyNucleus->GetOuterRadius()+projectile->GetOuterRadius();
      G4double aX=(2.*G4UniformRand()-1.)*impactMax;
      G4double aY=(2.*G4UniformRand()-1.)*impactMax;
      G4ThreeVector pos(aX, aY, -2.*impactMax);
      debug.push_back("Impact parameter");
      debug.push_back(aX);
      debug.push_back(aY);
      debug.push_back(-2.*impactMax);
      debug.dump();

      G4KineticTrackVector * initalState = new G4KineticTrackVector;
      G4LorentzVector cache = mom;
      cache.setZ(-cache.getZ());
      projectile->DoLorentzBoost(cache);
      projectile->StartLoop();
      G4Nucleon * aNuc;
      debug.push_back("Constituent energies after push");
      debug.push_back(mom);
      G4LorentzVector tmpV(0,0,0,0);
      G4LorentzVector nucleonMom(1./a1*mom);
      nucleonMom.setZ(nucleonMom.vect().mag());
      nucleonMom.setX(0);
      nucleonMom.setY(0);
      while( (aNuc=projectile->GetNextNucleon()) )
      {
	G4LorentzVector p4 = aNuc->GetMomentum();	
	tmpV+=p4;
	G4KineticTrack * it = new G4KineticTrack(aNuc, aNuc->GetPosition()+pos, nucleonMom );
	initalState->push_back(it);
      }
      debug.push_back(tmpV);
      debug.dump();
      result=theModel.Propagate(initalState, fancyNucleus);
      debug.push_back("################# Result size");
      debug.push_back(result->size());
      debug.dump();
      delete projectile;
      for_each(initalState->begin(), initalState->end(), Delete<G4KineticTrack>());
      delete initalState;
      
      if(result->size()==0) 
      {
        delete result; result=0;
        delete fancyNucleus;
      } 
      else
      {
        break;
      }     
    }
    debug.push_back("################# Through the loop ? "); debug.dump();
    
    //inverse transformation in case we swapped.
    fancyNucleus->StartLoop();  
    G4int resA(0), resZ(0); 
    G4Nucleon * aNuc;
    debug.push_back("getting at the hits"); debug.dump();
    while( (aNuc=fancyNucleus->GetNextNucleon()) )
    {
      debug.push_back("getting the hits"); debug.dump();
      if(!aNuc->AreYouHit())
      {
        resA++;
	resZ+=aNuc->GetDefinition()->GetPDGCharge();
      }
      else
      {
        debug.push_back(" ##### a hit ##### "); debug.dump();
      }
      debug.push_back("collected a hit"); debug.dump();
    }
    delete fancyNucleus;
    debug.push_back("have the hits"); 
    debug.push_back(resA);
    debug.push_back(resZ);
    debug.dump();
    // Calculate excitation energy
    G4LorentzVector iState = mom;
    iState.setT(iState.getT()+m2);

    G4LorentzVector fState(0,0,0,0);
    G4int i(0);
    for(i=0; i<result->size(); i++)
    {
      if( (*result)[i]->GetNewlyAdded() ) 
      {
        fState += G4LorentzVector( (*result)[i]->GetMomentum(), (*result)[i]->GetTotalEnergy() );
      }
    }
    G4LorentzVector momentum(iState-fState);

    //Make the fragment
    G4Fragment aProRes;
    aProRes.SetA(resA);
    aProRes.SetZ(resZ);
    aProRes.SetNumberOfParticles(0);
    aProRes.SetNumberOfCharged(0);
    aProRes.SetNumberOfHoles(a2-resA);
    aProRes.SetMomentum(momentum);

    // call precompound model
    G4ReactionProductVector * proFrag(0);
    proFrag = theProjectileFragmentation.DeExcite(aProRes);

    // collect the evaporation part
    G4ReactionProductVector::iterator ii;
    for(ii=proFrag->begin(); ii!=proFrag->end(); ii++)
    {
      result->push_back(*ii);
    }

    debug.push_back("################# done with evaporation"); debug.dump();
    // Rotate to lab
    G4LorentzRotation toZ;
    toZ.rotateZ(-1*mom.phi());
    toZ.rotateY(-1*mom.theta());
    G4LorentzRotation toLab(toZ.inverse());
  
    // Fill the particle change, while rotating. Boost from projectile breit-frame in case we swapped.  
    theParticleChange.Clear();
    theParticleChange.Initialize(aTrack);
    theParticleChange.SetStatusChange(fStopAndKill);
    theParticleChange.SetNumberOfSecondaries(result->size());
    for(i=0; i<result->size(); i++)
    {
      if((*result)[i]->GetNewlyAdded())
      {
	G4DynamicParticle * aNew = 
	new G4DynamicParticle((*result)[i]->GetDefinition(),
                              (*result)[i]->GetTotalEnergy(),
			      (*result)[i]->GetMomentum() );
	G4LorentzVector tmp = aNew->Get4Momentum();
	if(swapped)
	{
          tmp*=toBreit.inverse();
	}    
	tmp *= toLab;
	aNew->Set4Momentum(tmp);
	theParticleChange.AddSecondary(aNew);
      }
    }
    return &theResult;
  }
