//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//

#ifndef G4DecayStrongResonances_h
#define G4DecayStrongResonances_h 1

#include "G4Fancy3DNucleus.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4KineticTrackVector.hh"
#include "G4FragmentVector.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticleVector.hh"

#include <algorithm>

class G4DecayStrongResonances
{
public:
   G4DecayStrongResonances(){}      
   ~G4DecayStrongResonances(){}

private:
   G4int operator==(G4DecayStrongResonances& right) {return (this == &right);}
   G4int operator!=(G4DecayStrongResonances& right) {return (this != &right);}
   
   G4double theEnergy;
      
public:
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus)
   {
     // decay the strong resonances
     G4ReactionProductVector * theResult = new G4ReactionProductVector;
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
     
     size_t aResult=0;
     for (aResult=0; aResult < result1->size(); aResult++)
     {
       G4ParticleDefinition * pdef;
       pdef=result1->operator[](aResult)->GetDefinition();
       secondaries=NULL;
       if ( pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s )
       {
	  secondaries = result1->operator[](aResult)->Decay();
       }
       if ( secondaries == NULL )
       {
	  result->push_back(result1->operator[](aResult));
	  result1->operator[](aResult)=NULL;	//protect for clearAndDestroy 
       } 
       else
       {
	 for (size_t aSecondary=0; aSecondary<secondaries->size(); aSecondary++)
	 {
	     result1->push_back(secondaries->operator[](aSecondary));
	 }
	 delete secondaries;
       }
     }
     G4std::for_each(result1->begin(), result1->end(), DeleteKineticTrack());
     delete result1;
     
     // translate to ReactionProducts
     G4ReactionProduct * it = NULL;
     for(aResult=0; aResult < result->size(); aResult++)
     {
       it = new G4ReactionProduct();
       it->SetDefinition((*result)[aResult]->GetDefinition());
       it->SetMass((*result)[aResult]->GetDefinition()->GetPDGMass());
       it->SetTotalEnergy((*result)[aResult]->Get4Momentum().t());
       it->SetMomentum((*result)[aResult]->Get4Momentum().vect());
       
       theResult->push_back(it);
     }
     G4std::for_each(result->begin(), result->end(), DeleteKineticTrack());
     delete result;
     return theResult;
   }


private:   
};

#endif // G4DecayStrongResonances_h


