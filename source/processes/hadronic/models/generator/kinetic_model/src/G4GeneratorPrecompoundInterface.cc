#include "G4GeneratorPrecompoundInterface.hh"
#include "G4DynamicParticleVector.hh"
#include "G4IonTable.hh"

//
// HPW, 10DEC 98, the decay part originally written by Gunter Folger in his FTF-test-program.
//
      
   G4VParticleChange* G4GeneratorPrecompoundInterface::
   ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
   {
     cout << "G4GeneratorPrecompoundInterface: ApplyYourself interface called stand-allone."<< endl;
     cout << "This class is only a mediator between generator and precompound"<<endl;
     cout << "Please remove from your physics list."<<endl;
     G4Exception("SEVERE: G4GeneratorPrecompoundInterface model interface called stand-allone.");
     return new G4ParticleChange;
   }
   
   G4DynamicParticleVector* G4GeneratorPrecompoundInterface::
   Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
   {
     G4DynamicParticleVector * theTotalResult = new G4DynamicParticleVector;

     // decay the strong resonances
     G4KineticTrackVector *result1, *secondaries, *result;
     result1=theSecondaries;
     result=new G4KineticTrackVector();
     
     for (G4int aResult=0; aResult < result1->entries(); aResult++)
     {
       G4ParticleDefinition * pdef;
       pdef=result1->at(aResult)->GetDefinition();
       secondaries=NULL;
       if ( pdef->GetPDGWidth() > 0 || pdef->GetPDGLifeTime() < 1*ns )
       {
	  secondaries = result1->at(aResult)->Decay();
       }
       if ( secondaries == NULL )
       {
	  result->insert(result1->at(aResult));
	  result1->at(aResult)=NULL;	//protect for clearAndDestroy 
       } 
       else
       {
	 for (G4int aSecondary=0; aSecondary<secondaries->entries(); aSecondary++)
	 {
	     result1->append(secondaries->at(aSecondary));
	 }
	 delete secondaries;
       }
     }
     result1->clearAndDestroy();
     delete result1;
     
     
     // prepare the fragment
     G4Fragment anInitialState;
     G4int anA=theNucleus->GetMassNumber();
     G4int aZ=theNucleus->GetCharge();
     G4int numberOfEx = 0;
     G4int numberOfCh = 0;
     G4int numberOfHoles = 0;
     G4double exEnergy = 0;
     G4ThreeVector exciton3Momentum(0,0,0);
     // loop over secondaries
     for(G4int list=0; list < result->entries(); list++)
     {
       G4KineticTrack *aTrack = result->at(list);
       if(aTrack->GetDefinition() != G4Proton::Proton() && 
          aTrack->GetDefinition() != G4Neutron::Neutron())
       {
         theTotalResult->insert(new G4DynamicParticle(aTrack->GetDefinition(), aTrack->Get4Momentum()));            
       }
       else if(aTrack->Get4Momentum().t() - aTrack->Get4Momentum().mag()>80*MeV)
       {
         theTotalResult->insert(new G4DynamicParticle(aTrack->GetDefinition(), aTrack->Get4Momentum()));            
       }
       else if(aTrack->GetPosition().mag() > theNucleus->GetNuclearRadius())
       {
         theTotalResult->insert(new G4DynamicParticle(aTrack->GetDefinition(), aTrack->Get4Momentum()));            
       }
       else
       {
         // within the nucleus, neutron or proton
         // now calculate  A, Z of the fragment, momentum, number of exciton states
         anA++;;
         numberOfEx++;
         aZ += aTrack->GetDefinition()->GetPDGCharge();
         numberOfCh += aTrack->GetDefinition()->GetPDGCharge();
         exciton3Momentum += aTrack->Get4Momentum().vect();
         exEnergy += (aTrack->Get4Momentum().t()-aTrack->Get4Momentum().m());
       }
     }
     
     // loop over wounded nucleus
     G4Nucleon * theCurrentNucleon = theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : NULL;
     while(theCurrentNucleon != NULL)
     {
       if(theCurrentNucleon->AreYouHit()) 
       {
         numberOfHoles++;
         numberOfEx++;
         anA--;
         aZ -= theCurrentNucleon->GetDefinition()->GetPDGCharge();
         exciton3Momentum -= theCurrentNucleon->Get4Momentum().vect();
         exEnergy+=theCurrentNucleon->GetBindingEnergy();
       }
       theCurrentNucleon = theNucleus->GetNextNucleon();
     }   
     
     G4double residualMass =  
              G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ ,anA);
     residualMass += exEnergy;
     G4LorentzVector exciton4Momentum(exciton3Momentum, 
                                      sqrt(exciton3Momentum.mag2()+residualMass*residualMass));
     
     anInitialState.SetA(anA);
     anInitialState.SetZ(aZ);
     anInitialState.SetNumberOfCharged(numberOfCh);
     anInitialState.SetNumberOfHoles(numberOfHoles);
     anInitialState.SetNumberOfExcitons(numberOfEx);
     anInitialState.SetMomentum(exciton4Momentum);
     anInitialState.SetExcitationEnergy(exEnergy);

     // call pre-compound
     const G4Fragment aFragment(anInitialState);
     G4DynamicParticleVector * aPreResult = theDeExcitation->DeExcite(aFragment);
   
     // fill pre-compound part into the result, and return
     for(G4int ll=0; ll<aPreResult->entries(); ll++)
     {
       theTotalResult->insert(aPreResult->at(ll));
     }
     delete aPreResult;
     
     // now return
     return theTotalResult;
   }
