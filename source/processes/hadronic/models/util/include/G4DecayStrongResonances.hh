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
#include "G4HadFinalState.hh"
#include "G4DynamicParticleVector.hh"
#include "G4HadTmpUtil.hh"

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
   G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* )
   {
     // decay the strong resonances
     //static int call_count = 0;
     //if(call_count++<10)
     //{
     //  G4cout << "Security print-out: Entering G4DecayStrongResonances::Propagate";
     //}
     G4ReactionProductVector * theResult;
     try
     {
       theResult = new G4ReactionProductVector;
     }
     catch(...)
     {
       throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: out of memory ");
     }
     G4KineticTrackVector *result1, *secondaries, *result;
     if(!theSecondaries)
     {
       throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: 0x0 input vector ");
     }
     result1=theSecondaries;
     try
     {
       result=new G4KineticTrackVector();
     }
     catch(...)
     {
       throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: out of memory in ");
     }
          
     size_t aResult=0;
     for (aResult=0; aResult < result1->size(); aResult++)
     {
       G4ParticleDefinition * pdef;
       if(!result1->operator[](aResult))
       {
         throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: null pointer in input vector!!! ");
       }
       pdef=result1->operator[](aResult)->GetDefinition();
       if(!pdef)
       {
        throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: 0x0 particle definition ");
       }
       secondaries=NULL;
       if ( pdef->GetPDGWidth() > 0 && pdef->GetPDGLifeTime() < 5E-17*s )
       {
	  try
	  {
 	    secondaries = result1->operator[](aResult)->Decay();
	  }
          catch(...)
          {
            throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: failing in Decay ");
          }
       }
       if ( secondaries == NULL )
       {
	  try
	  {
	    result->push_back(result1->operator[](aResult));
	  }
          catch(...)
          {
            throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: push_back failed - out of memory ");
          }
	  result1->operator[](aResult)=NULL;	//protect for clearAndDestroy 
       } 
       else
       {
	 for (size_t aSecondary=0; aSecondary<secondaries->size(); aSecondary++)
	 {
	   try
	   {
	     result1->push_back(secondaries->operator[](aSecondary));
	   }
           catch(...)
           {
             throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: push_back  1 failed - out of mem");
           }
	 }
	 if(secondaries) delete secondaries;
       }
     }
     try
     {
       std::for_each(result1->begin(), result1->end(), G4Delete());
       delete result1;
     }
     catch(...)
     {
       throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: memory corruption.");
     }
     
     // translate to ReactionProducts
     G4ReactionProduct * it = NULL;
     for(aResult=0; aResult < result->size(); aResult++)
     {
       try
       {
         it = new G4ReactionProduct();
       }
       catch(...)
       {
          throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: out of memory ");
       }
       it->SetDefinition((*result)[aResult]->GetDefinition());
       it->SetMass((*result)[aResult]->GetDefinition()->GetPDGMass());
       it->SetTotalEnergy((*result)[aResult]->Get4Momentum().t());
       it->SetMomentum((*result)[aResult]->Get4Momentum().vect());
       
       try
       {
         theResult->push_back(it);
       }
       catch(...)
       {
          throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: push to result failed - out of mem.");
       }
     }
     try
     {
       std::for_each(result->begin(), result->end(), G4Delete());
       delete result;
     }
     catch(...)
     {
       throw G4HadronicException(__FILE__, __LINE__, "DecayStrongRes: memory corruption at end.");
     }
     return theResult;
   }
};

#endif // G4DecayStrongResonances_h


