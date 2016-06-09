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
#ifndef G4HadronQEDBuilder_h
#define G4HadronQEDBuilder_h 1

#include "globals.hh"
#include "G4ios.hh"

#include "G4MultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ProcessManager.hh"

#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4SigmaMinus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4SigmaPlus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4XiMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4OmegaMinus.hh" 
#include "G4AntiOmegaMinus.hh"
#include "plist.tmp"
#include "ParticleCodeMap.hh"

class G4HadronQEDBuilder 
{
  public: 
    G4HadronQEDBuilder();
    virtual ~G4HadronQEDBuilder();

  public: 
    virtual void Build();

  private:  
    void RegisterOne(G4ProcessManager* aP, G4MultipleScattering * aM, G4hIonisation* aI);

   public: // to get gcc 2.95 satisfied...
   typedef Plist<G4PionPlus, G4PionMinus, G4KaonPlus, G4KaonMinus, G4Proton, G4AntiProton, G4SigmaMinus,
          G4AntiSigmaMinus, G4SigmaPlus, G4AntiSigmaPlus, G4XiMinus, G4AntiXiMinus, G4OmegaMinus,
	  G4AntiOmegaMinus> theParticles;
	  
   typedef std::pair<G4ParticleDefinition *, std::pair<G4MultipleScattering *, G4hIonisation *> > entryType;
   
   static std::vector< entryType > & theCache()
   {
     static std::vector< entryType > cache;
     return cache;
   }
   private:

   G4bool wasActivated;
   
   struct Clear
   { void operator () (entryType & aE)
     {
       G4ProcessManager * aP = aE.first->GetProcessManager();
       if(aP) aP->RemoveProcess(aE.second.first);
       if(aP) aP->RemoveProcess(aE.second.second);
     }
   };
   
   struct Register
   {
     template <class T> 
     struct Fun // to get gcc 2.95 satisfied...
     {
       void operator()()
       {
	 G4ParticleDefinition * aPart = 
              G4ParticleTable::GetParticleTable()->FindParticle(ParticleCodeMap<T>::PDG_CODE);
	 G4cout << aPart<<G4endl;
	 // G4cout << "===== > Calling for "<<aPart->GetParticleName()<<G4endl;
	 G4ProcessManager * aP = aPart->GetProcessManager();
	 G4hIonisation * aI = new G4hIonisation;
	 G4MultipleScattering * aM = new G4MultipleScattering;
	 aP->AddProcess(aI, ordInActive,2, 2);
	 aP->AddProcess(aM);
	 aP->SetProcessOrdering(aM, idxAlongStep, 1);
	 aP->SetProcessOrdering(aM, idxPostStep, 1);
	 G4HadronQEDBuilder::theCache().push_back(std::make_pair(aPart, std::pair<G4MultipleScattering *, G4hIonisation *>(aM, aI) ) );
       }
     };
   };
};

// 2002 by J.P. Wellisch

#endif

