
#ifndef G4StringInfoDump_h
#define G4StringInfoDump_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4VIntraNuclearTransportModel.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticleVector.hh"

class G4StringInfoDump : public G4VIntraNuclearTransportModel 
{
public:
   G4StringInfoDump(){}      
   G4StringInfoDump(G4double anEnergy)
   {
     theEnergy = anEnergy;
   }      
   ~G4StringInfoDump(){}

private:
   G4int operator==(G4StringInfoDump& right) {return (this == &right);}
   G4int operator!=(G4StringInfoDump& right) {return (this != &right);}
   
   G4double theEnergy;
      
public:
   G4VParticleChange* ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
   { return new G4ParticleChange;}

   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
   {
     // decay the strong resonances
     G4KineticTrackVector *result1, *secondaries, *result;
     result1=theSecondaries;
     result=new G4KineticTrackVector();
     
     G4ThreeVector theCurrentVelocity;
     theCurrentVelocity.setX(0);    
     theCurrentVelocity.setY(0); 
     G4double vz = theEnergy/(sqrt(theEnergy*theEnergy + 
                     G4Proton::Proton()->GetPDGMass()*G4Proton::Proton()->GetPDGMass()) 
                   + G4Proton::Proton()->GetPDGMass());
     theCurrentVelocity.setZ(vz);

     G4double preDecayEtWithinCuts = 0;
     double totalEnergy = 0;
     int pos, neg, ch, al;
     pos = neg = ch = 0;
     al = result1->entries();
     
     // some info about the energy of particles within the nucleus
     G4double nuclearRadius = theNucleus->GetNuclearRadius(25*perCent);
     G4double insideEnergy = 0;
     for (G4int aRes=0; aRes < result1->entries(); aRes++)
     {
       G4cout << "PREDECAYDUMP ";
       G4cout << result1->at(aRes)->GetDefinition()->GetPDGCharge()<<" ";
       G4cout << result1->at(aRes)->GetDefinition()->GetPDGEncoding()<<" ";
       G4cout << result1->at(aRes)->Get4Momentum()<<" ";
       G4cout << result1->at(aRes)->GetPosition()<<" ";
       totalEnergy += result1->at(aRes)->Get4Momentum().t();
       if(result1->at(aRes)->GetPosition().mag()<nuclearRadius)
       {
         insideEnergy +=result1->at(aRes)->Get4Momentum().t();
       }
       if(result1->at(aRes)->GetDefinition()->GetPDGCharge()>0)
       {
         pos++;
         ch++;
       }
       if(result1->at(aRes)->GetDefinition()->GetPDGCharge()<0)
       {
         neg++;
         ch++;
       }
       G4LorentzVector Mom = result1->at(aRes)->Get4Momentum();
       double currentEta = Mom.rapidity();
       if(currentEta>-0.1&&currentEta<2.9) 
	      preDecayEtWithinCuts += Mom.perp();
       Mom.boost(-theCurrentVelocity);
       G4cout << Mom << G4endl;
     }
     G4cout << "PREDECAYETWITHINCUTS "<<preDecayEtWithinCuts<<" "
            <<pos<<" "<<neg<<" "<<ch<<" "<<al
	    << " total energy = "<< totalEnergy 
	    << " inside energy = "<<insideEnergy << G4endl;
     totalEnergy = 0;
     for (G4int aResult=0; aResult < result1->entries(); aResult++)
     {
       G4ParticleDefinition * pdef;
       pdef=result1->at(aResult)->GetDefinition();
       secondaries=NULL;
       if ( pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s )
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
     double etWithinCuts = 0;
     int positives, negatives, charged, all;
     positives = negatives = charged = 0;
     all = result->length();
     for(G4int i=0; i<result->length(); i++)
     {
       totalEnergy += result->at(i)->Get4Momentum().t();
       G4cout << "STRINGDUMP ";
       G4cout << result->at(i)->GetDefinition()->GetPDGCharge()<<" ";
       G4cout << result->at(i)->GetDefinition()->GetPDGEncoding()<<" ";
       G4cout << result->at(i)->Get4Momentum()<<" ";
       G4cout << result->at(i)->GetPosition()<<" ";
       if(abs(result->at(i)->GetDefinition()->GetPDGCharge())>1.1)
	 G4Exception("multiply charged particle after decay");
      if(result->at(i)->GetDefinition()->GetPDGCharge()>0.5)
       {
         positives++;
         charged++;
       }
       if(result->at(i)->GetDefinition()->GetPDGCharge()<-0.5)
       {
         negatives++;
         charged++;
       }
       if(abs(result->at(i)->GetDefinition()->GetPDGCharge())>1.5) 
           G4Exception("still multiply charged after decay");
       
       G4LorentzVector Mom = result->at(i)->Get4Momentum();
       double currentEta = Mom.rapidity();
       if(currentEta>-0.1&&currentEta<2.9) 
       {
	 const G4ParticleDefinition * theCurrentDef = result->at(i)->GetDefinition();
	 if(theCurrentDef->GetBaryonNumber()>0)
	 {
	   double energy = Mom.t() - G4Proton::Proton()->GetPDGMass();
	   etWithinCuts += sin(Mom.theta())*energy;
	 }
	 else if(theCurrentDef->GetBaryonNumber()<0)
	 {
	   double energy = Mom.t() + G4Proton::Proton()->GetPDGMass();
	   etWithinCuts += sin(Mom.theta())*energy;
	 }
	 else
	 {
	   double energy = Mom.t();
	   etWithinCuts += sin(Mom.theta())*energy;
	 }
       }
       Mom.boost(-theCurrentVelocity);
       G4cout << Mom << G4endl;
     }
     G4cout << "ETWITHINCUTS "<<etWithinCuts<<" "<<positives<<" "
     <<negatives<<" "<<charged<<" "<<all<< " total energy = "<< totalEnergy <<G4endl;
     G4ReactionProductVector * theFinalResult = new G4ReactionProductVector;
     G4ReactionProduct * theSec;
     for(G4int it=0; it<result->length(); it++)
     {
       theSec = new G4ReactionProduct(result->at(it)->GetDefinition());
       G4LorentzVector current4Mom = result->at(it)->Get4Momentum();
       theSec->SetTotalEnergy(current4Mom.t());
       theSec->SetMomentum(current4Mom.vect());
       theFinalResult->insert(theSec);
     }
     result->clearAndDestroy();
     delete result;
     return theFinalResult;
   }


private:   
};

#endif // G4StringInfoDump_h


