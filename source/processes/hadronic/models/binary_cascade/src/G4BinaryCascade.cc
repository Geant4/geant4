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

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4LorentzRotation.hh"
#include "G4BinaryCascade.hh"
#include "G4KineticTrackVector.hh"
#include "G4DecayKineticTracks.hh"
#include "G4ReactionProductVector.hh"
#include "G4Track.hh"
#include "G4V3DNucleus.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4Scatterer.hh"
#include "G4MesonAbsorption.hh"
#include "G4ping.hh"
#include "G4Delete.hh"

#include "G4CollisionManager.hh"
#include "G4Absorber.hh"

#include "G4CollisionInitialState.hh"
#include "G4ListOfCollisions.hh"
#include "G4Fragment.hh"
#include "G4RKPropagation.hh"

#include "G4NuclearShellModelDensity.hh"
#include "G4NuclearFermiDensity.hh"
#include "G4FermiMomentum.hh"

#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4HadronicInteractionRegistry.hh"

#include "G4FermiPhaseSpaceDecay.hh"

#include "G4PreCompoundModel.hh"

#include <algorithm>
#include "G4ShortLivedConstructor.hh"
#include <typeinfo>

//   turn on general debugging info, and consistency checks

//#define debug_G4BinaryCascade 1

//  more detailed debugging -- deprecated
//#define debug_H1_BinaryCascade 1

//  specific debugging info per method or functionality
//#define debug_BIC_ApplyCollision 1
//#define debug_BIC_CheckPauli 1
//#define debug_BIC_CorrectFinalPandE 1
//#define debug_BIC_Propagate 1
//#define debug_BIC_Propagate_Excitation 1
//#define debug_BIC_Propagate_Collisions 1
//#define debug_BIC_Propagate_finals 1
//#define debug_BIC_DoTimeStep 1
//#define debug_BIC_CorrectBarionsOnBoundary 1
//#define debug_BIC_GetExcitationEnergy 1
//#define debug_BIC_DeexcitationProducts 1
//#define debug_BIC_FinalNucleusMomentum 1
//#define debug_BIC_Final4Momentum 1
//#define debug_BIC_FillVoidnucleus 1
//#define debug_BIC_FindFragments 1
//#define debug_BIC_BuildTargetList 1
//#define debug_BIC_FindCollision 1
//#define debug_BIC_return 1

//-------
//#if defined(debug_G4BinaryCascade)
#if 0
  #define _CheckChargeAndBaryonNumber_(val) CheckChargeAndBaryonNumber(val)
  //#define debugCheckChargeAndBaryonNumberverbose 1
#else
  #define _CheckChargeAndBaryonNumber_(val)
#endif
//#if defined(debug_G4BinaryCascade)
#if 0
  #define _DebugEpConservation(val)  DebugEpConservation(val)
  //#define debugCheckChargeAndBaryonNumberverbose 1
#else
  #define _DebugEpConservation(val)
#endif

#ifdef G4MULTITHREADED
   G4Mutex G4BinaryCascade::BICMutex = G4MUTEX_INITIALIZER;
#endif
   G4int G4BinaryCascade::theBIC_ID = -1;

//
//  C O N S T R U C T O R S   A N D   D E S T R U C T O R S
//
G4BinaryCascade::G4BinaryCascade(G4VPreCompoundModel* ptr) :
G4VIntraNuclearTransportModel("Binary Cascade", ptr)
{
    // initialise the resonance sector
    G4ShortLivedConstructor ShortLived;
    ShortLived.ConstructParticle();

    theCollisionMgr = new G4CollisionManager;
    theDecay=new G4BCDecay;
    theImR.push_back(theDecay);
    theLateParticle= new G4BCLateParticle;
    G4MesonAbsorption * aAb=new G4MesonAbsorption;
    theImR.push_back(aAb);
    G4Scatterer * aSc=new G4Scatterer;
    theH1Scatterer = new G4Scatterer;
    theImR.push_back(aSc);

    thePropagator = new G4RKPropagation;
    theCurrentTime = 0.;
    theBCminP = 45*MeV;
    theCutOnP = 90*MeV;
    theCutOnPAbsorb= 0*MeV;   // No Absorption of slow Mesons, other than above G4MesonAbsorption

    // reuse existing pre-compound model
    if(!ptr) {
      G4HadronicInteraction* p =
	G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
      G4VPreCompoundModel* pre = static_cast<G4VPreCompoundModel*>(p);
      if(!pre) { pre = new G4PreCompoundModel(); }
      SetDeExcitation(pre);
    }
    theExcitationHandler = GetDeExcitation()->GetExcitationHandler();
    SetMinEnergy(0.0*GeV);
    SetMaxEnergy(10.1*GeV);
    //PrintWelcomeMessage();
    thePrimaryEscape = true;
    thePrimaryType = 0;

    SetEnergyMomentumCheckLevels(1.0*perCent, 1.0*MeV);

    // init data members
    currentA=currentZ=0;
    lateA=lateZ=0;
    initialA=initialZ=0;
    projectileA=projectileZ=0;
    currentInitialEnergy=initial_nuclear_mass=0.;
    massInNucleus=0.;
    theOuterRadius=0.;
    if ( theBIC_ID == -1 ) {
#ifdef G4MULTITHREADED
       G4MUTEXLOCK(&G4BinaryCascade::BICMutex);
       if ( theBIC_ID == -1 ) {
#endif
    	   theBIC_ID = G4PhysicsModelCatalog::Register("Binary Cascade");
#ifdef G4MULTITHREADED
       }
       G4MUTEXUNLOCK(&G4BinaryCascade::BICMutex);
#endif
    }

}

/*
G4BinaryCascade::G4BinaryCascade(const G4BinaryCascade& )
: G4VIntraNuclearTransportModel("Binary Cascade")
{
}
 */

G4BinaryCascade::~G4BinaryCascade()
{
    ClearAndDestroy(&theTargetList);
    ClearAndDestroy(&theSecondaryList);
    ClearAndDestroy(&theCapturedList);
    delete thePropagator;
    delete theCollisionMgr;
    std::for_each(theImR.begin(), theImR.end(), Delete<G4BCAction>());
    delete theLateParticle;
    //delete theExcitationHandler;
    delete theH1Scatterer;
}

void G4BinaryCascade::ModelDescription(std::ostream& outFile) const
{
    outFile << "G4BinaryCascade is an intra-nuclear cascade model in which\n"
            << "an incident hadron collides with a nucleon, forming two\n"
            << "final-state particles, one or both of which may be resonances.\n"
            << "The resonances then decay hadronically and the decay products\n"
            << "are then propagated through the nuclear potential along curved\n"
            << "trajectories until they re-interact or leave the nucleus.\n"
            << "This model is valid for incident pions up to 1.5 GeV and\n"
            << "nucleons up to 10 GeV.\n"
            << "The remaining excited nucleus is handed on to ";
            if (theDeExcitation)                // pre-compound
            {
              outFile << theDeExcitation->GetModelName() << " : \n ";
              theDeExcitation->DeExciteModelDescription(outFile);
            }
            else if (theExcitationHandler)    // de-excitation
            {
               outFile << "G4ExcitationHandler";    //theExcitationHandler->GetModelName();
               theExcitationHandler->ModelDescription(outFile);
            }
            else
            {
               outFile << "void.\n";
            }
    outFile<< " \n";
}
void G4BinaryCascade::PropagateModelDescription(std::ostream& outFile) const
{
    outFile << "G4BinaryCascade propagtes secondaries produced by a high\n"
            << "energy model through the wounded nucleus.\n"
            << "Secondaries are followed after the formation time and if\n"
            << "within the nucleus are propagated through the nuclear\n"
            << "potential along curved trajectories until they interact\n"
            << "with a nucleon, decay, or leave the nucleus.\n"
            << "An interaction of a secondary with a nucleon produces two\n"
            << "final-state particles, one or both of which may be resonances.\n"
            << "Resonances decay hadronically and the decay products\n"
            << "are in turn propagated through the nuclear potential along curved\n"
            << "trajectories until they re-interact or leave the nucleus.\n"
            << "This model is valid for pions up to 1.5 GeV and\n"
            << "nucleons up to about 3.5 GeV.\n"
            << "The remaining excited nucleus is handed on to ";
    if (theDeExcitation)                // pre-compound
    {
      outFile << theDeExcitation->GetModelName() << " : \n ";
      theDeExcitation->DeExciteModelDescription(outFile);
    }
    else if (theExcitationHandler)    // de-excitation
    {
       outFile << "G4ExcitationHandler";    //theExcitationHandler->GetModelName();
       theExcitationHandler->ModelDescription(outFile);
    }
    else
    {
       outFile << "void.\n";
    }
outFile<< " \n";
}

//----------------------------------------------------------------------------

//
//      I M P L E M E N T A T I O N
//


//----------------------------------------------------------------------------
G4HadFinalState * G4BinaryCascade::ApplyYourself(const G4HadProjectile & aTrack,
        G4Nucleus & aNucleus)
//----------------------------------------------------------------------------
{
    if(getenv("BCDEBUG") ) G4cerr << " ######### Binary Cascade Reaction starts ######### "<< G4endl;

    G4LorentzVector initial4Momentum = aTrack.Get4Momentum();
    const G4ParticleDefinition * definition = aTrack.GetDefinition();

    if(initial4Momentum.e()-initial4Momentum.m()<theBCminP &&
            ( definition==G4Neutron::NeutronDefinition() || definition==G4Proton::ProtonDefinition() ) )
    {
        return theDeExcitation->ApplyYourself(aTrack, aNucleus);
    }

    theParticleChange.Clear();
    // initialize the G4V3DNucleus from G4Nucleus
    the3DNucleus = new G4Fancy3DNucleus;

    // Build a KineticTrackVector with the G4Track
    G4KineticTrackVector * secondaries;// = new G4KineticTrackVector;
    G4ThreeVector initialPosition(0., 0., 0.); // will be set later

    if(!getenv("I_Am_G4BinaryCascade_Developer") )
    {
        if(definition!=G4Neutron::NeutronDefinition() &&
                definition!=G4Proton::ProtonDefinition() &&
                definition!=G4PionPlus::PionPlusDefinition() &&
                definition!=G4PionMinus::PionMinusDefinition() )
        {
            G4cerr << "You are trying to use G4BinaryCascade with " <<definition->GetParticleName()<<" as projectile."<<G4endl;
            G4cerr << "G4BinaryCascade should not be used for projectiles other than nucleons or pions."<<G4endl;
            G4cerr << "If you want to continue, please switch on the developer environment: "<<G4endl;
            G4cerr << "setenv I_Am_G4BinaryCascade_Developer 1 "<<G4endl<<G4endl;
            throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade - used for unvalid particle type - Fatal");
        }
    }

    // keep primary
    thePrimaryType = definition;
    thePrimaryEscape = false;

    G4double timePrimary=aTrack.GetGlobalTime();

    // try until an interaction will happen
    G4ReactionProductVector * products=0;
    G4int interactionCounter = 0,collisionLoopMaxCount;
    do
    {
        // reset status that could be changed in previous loop event
        theCollisionMgr->ClearAndDestroy();

        if(products != 0)
        {  // free memory from previous loop event
            ClearAndDestroy(products);
            delete products;
            products=0;
        }

        G4int massNumber=aNucleus.GetA_asInt();
        the3DNucleus->Init(massNumber, aNucleus.GetZ_asInt());
        thePropagator->Init(the3DNucleus);
        G4KineticTrack * kt;
		  collisionLoopMaxCount = 200;
        do                  // sample impact parameter until collisions are found
        {
            theCurrentTime=0;
            G4double radius = the3DNucleus->GetOuterRadius()+3*fermi;
            initialPosition=GetSpherePoint(1.1*radius, initial4Momentum);  // get random position
            kt = new G4KineticTrack(definition, 0., initialPosition, initial4Momentum);
            kt->SetState(G4KineticTrack::outside);
            // secondaries has been cleared by Propagate() in the previous loop event
            secondaries= new G4KineticTrackVector;
            secondaries->push_back(kt);
            if(massNumber > 1) // 1H1 is special case
            {
                products = Propagate(secondaries, the3DNucleus);
            } else {
                products = Propagate1H1(secondaries,the3DNucleus);
            }
				    // until we FIND a collision ... or give up
        } while(! products && --collisionLoopMaxCount>0);   /* Loop checking, 31.08.2015, G.Folger */

        if(++interactionCounter>99) break;
		    // ...until we find an ALLOWED collision ... or give up
    } while(products && products->size() == 0);   /* Loop checking, 31.08.2015, G.Folger */

    if(products && products->size()>0)
    {
        //  G4cout << "BIC Applyyourself: number of products " << products->size() << G4endl;

        // Fill the G4ParticleChange * with products
        theParticleChange.SetStatusChange(stopAndKill);
        G4ReactionProductVector::iterator iter;

        for(iter = products->begin(); iter != products->end(); ++iter)
        {
        	G4DynamicParticle * aNewDP =
                    new G4DynamicParticle((*iter)->GetDefinition(),
                            (*iter)->GetTotalEnergy(),
                            (*iter)->GetMomentum());
        	G4HadSecondary aNew = G4HadSecondary(aNewDP);
            G4double time=(*iter)->GetFormationTime();
            if(time < 0.0) { time = 0.0; }
            aNew.SetTime(timePrimary + time);
            aNew.SetCreatorModelType((*iter)->GetCreatorModel());
            theParticleChange.AddSecondary(aNew);
        }

         //DebugFinalEpConservation(aTrack, products);


    } else {  // no interaction, return primary
        if(getenv("BCDEBUG") ) G4cerr << " ######### Binary Cascade Reaction void, return intial state ######### "<< G4endl;
        theParticleChange.SetStatusChange(isAlive);
        theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
        theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    }

    if ( products ) 
	 {
	    ClearAndDestroy(products);
       delete products;
    }
	 
    delete the3DNucleus;
    the3DNucleus = NULL;

    if(getenv("BCDEBUG") ) G4cerr << " ######### Binary Cascade Reaction ends ######### "<< G4endl;

    return &theParticleChange;
}
//----------------------------------------------------------------------------
G4ReactionProductVector * G4BinaryCascade::Propagate(
        G4KineticTrackVector * secondaries, G4V3DNucleus * aNucleus)
//----------------------------------------------------------------------------
{
    G4ping debug("debug_G4BinaryCascade");
#ifdef debug_BIC_Propagate
    G4cout << "G4BinaryCascade Propagate starting -------------------------------------------------------" <<G4endl;
#endif

    the3DNucleus=aNucleus;
    G4ReactionProductVector * products = new G4ReactionProductVector;
    theOuterRadius = the3DNucleus->GetOuterRadius();
    theCurrentTime=0;
    theProjectile4Momentum=G4LorentzVector(0,0,0,0);
    theMomentumTransfer=G4ThreeVector(0,0,0);
    // build theSecondaryList, theProjectileList and theCapturedList
    ClearAndDestroy(&theCapturedList);
    ClearAndDestroy(&theSecondaryList);
    theSecondaryList.clear();
    ClearAndDestroy(&theFinalState);
    std::vector<G4KineticTrack *>::iterator iter;
    theCollisionMgr->ClearAndDestroy();

    theCutOnP=90*MeV;
    if(the3DNucleus->GetMass()>30) theCutOnP = 70*MeV;
    if(the3DNucleus->GetMass()>60) theCutOnP = 50*MeV;
    if(the3DNucleus->GetMass()>120) theCutOnP = 45*MeV;


    BuildTargetList();

#ifdef debug_BIC_GetExcitationEnergy
    G4cout << "ExcitationEnergy0 " << GetExcitationEnergy() << G4endl;
#endif

    thePropagator->Init(the3DNucleus);

    G4bool success = BuildLateParticleCollisions(secondaries);
    if (! success )   // fails if no excitation energy left....
    {
       products=HighEnergyModelFSProducts(products, secondaries);
       ClearAndDestroy(secondaries);
       delete secondaries;

#ifdef debug_G4BinaryCascade
       G4cout << "G4BinaryCascade::Propagate: warning - high energy model failed energy conservation, returning unchanged high energy final state" << G4endl;
#endif

       return products;
    }
    // check baryon and charge ...

    _CheckChargeAndBaryonNumber_("lateparticles");
    _DebugEpConservation(" be4 findcollisions");

    // if called stand alone find first collisions
    FindCollisions(&theSecondaryList);


    if(theCollisionMgr->Entries() == 0 )      //late particles ALWAYS create Entries
    {
        //G4cout << " no collsions -> return 0" << G4endl;
        delete products;
#ifdef debug_BIC_return
        G4cout << "return @ begin2,  no collisions "<< G4endl;
#endif
        return 0;
    }

    // end of initialization: do the job now
    // loop until there are no more collisions


    G4bool haveProducts = false;
    G4int collisionCount=0;
	 G4int collisionLoopMaxCount=1000000;
    while(theCollisionMgr->Entries() > 0 && currentZ && --collisionLoopMaxCount>0)  /* Loop checking, 31.08.2015, G.Folger */  
    {
      if(Absorb()) {  // absorb secondaries, pions only
            haveProducts = true;
      }
      if(Capture()) { // capture secondaries, nucleons only
            haveProducts = true;
        }

        // propagate to the next collision if any (collisions could have been deleted
        // by previous absorption or capture)
        if(theCollisionMgr->Entries() > 0)
        {
            G4CollisionInitialState *
            nextCollision = theCollisionMgr->GetNextCollision();
#ifdef debug_BIC_Propagate_Collisions
            G4cout << " NextCollision  * , Time, curtime = " << nextCollision << " "
                    <<nextCollision->GetCollisionTime()<< " " <<
                    theCurrentTime<< G4endl;
#endif
            if (!DoTimeStep(nextCollision->GetCollisionTime()-theCurrentTime) )
            {
                // Check if nextCollision is still valid, ie. particle did not leave nucleus
                if (theCollisionMgr->GetNextCollision() != nextCollision )
                {
                    nextCollision = 0;
                }
            }
           //_DebugEpConservation("Stepped");

            if( nextCollision )
            {
                if (ApplyCollision(nextCollision))
                {
                    //G4cerr << "ApplyCollision success " << G4endl;
                    haveProducts = true;
                    collisionCount++;
                    //_CheckChargeAndBaryonNumber_("ApplyCollision");
                    //_DebugEpConservation("ApplyCollision");

                } else {
                    //G4cerr << "ApplyCollision failure " << G4endl;
                    theCollisionMgr->RemoveCollision(nextCollision);
                }
            }
        }
    }

    //--------- end of on Collisions
    //G4cout << "currentZ @ end loop " << currentZ << G4endl;
    G4int nProtons(0);
    for(iter = theTargetList.begin(); iter != theTargetList.end(); ++iter)
    {
    	if ( (*iter)->GetDefinition() == G4Proton::Proton() ) ++nProtons;
    }
    if ( ! theTargetList.size() || ! nProtons ){
        // nucleus completely destroyed, fill in ReactionProductVector
       products = FillVoidNucleusProducts(products);
#ifdef debug_BIC_return
        G4cout << "return @ Z=0 after collision loop "<< G4endl;
        PrintKTVector(&theSecondaryList,std::string(" theSecondaryList"));
        G4cout << "theTargetList size: " << theTargetList.size() << G4endl;
        PrintKTVector(&theTargetList,std::string(" theTargetList"));
        PrintKTVector(&theCapturedList,std::string(" theCapturedList"));

        G4cout << " ExcitE be4 Correct : " <<GetExcitationEnergy() << G4endl;
        G4cout << " Mom Transfered to nucleus : " << theMomentumTransfer << " " << theMomentumTransfer.mag() << G4endl;
        PrintKTVector(&theFinalState,std::string(" FinalState uncorrected"));
        G4cout << "returned products: " << products->size() << G4endl;
        _CheckChargeAndBaryonNumber_("destroyed Nucleus");
        _DebugEpConservation("destroyed Nucleus");
#endif

        return products;
    }

    // No more collisions: absorb, capture and propagate the secondaries out of the nucleus
    if(Absorb()) {
        haveProducts = true;
        // G4cout << "Absorb sucess " << G4endl;
    }

    if(Capture()) {
        haveProducts = true;
        // G4cout << "Capture sucess " << G4endl;
    }

    if(!haveProducts)  // no collisions happened. Return an empty vector.
    {
#ifdef debug_BIC_return
        G4cout << "return 3, no products "<< G4endl;
#endif
        return products;
    }


#ifdef debug_BIC_Propagate
    G4cout << " Momentum transfer to Nucleus " << theMomentumTransfer << " " << theMomentumTransfer.mag() << G4endl;
    G4cout << "  Stepping particles out...... " << G4endl;
#endif

    StepParticlesOut();
    _DebugEpConservation("stepped out");


    if ( theSecondaryList.size() > 0 )
    {
#ifdef debug_G4BinaryCascade
        G4cerr << "G4BinaryCascade: Warning, have active particles at end" << G4endl;
        PrintKTVector(&theSecondaryList, "active particles @ end  added to theFinalState");
#endif
        //  add left secondaries to FinalSate
        for ( iter =theSecondaryList.begin(); iter != theSecondaryList.end(); ++iter)
        {
            theFinalState.push_back(*iter);
        }
        theSecondaryList.clear();

    }
    while ( theCollisionMgr->Entries() > 0 )                /* Loop checking, 31.08.2015, G.Folger */
    {
#ifdef debug_G4BinaryCascade
        G4cerr << " Warning: remove left over collision(s) " << G4endl;
#endif
        theCollisionMgr->RemoveCollision(theCollisionMgr->GetNextCollision());
    }

#ifdef debug_BIC_Propagate_Excitation

    PrintKTVector(&theSecondaryList,std::string(" theSecondaryList"));
    G4cout << "theTargetList size: " << theTargetList.size() << G4endl;
    //  PrintKTVector(&theTargetList,std::string(" theTargetList"));
    PrintKTVector(&theCapturedList,std::string(" theCapturedList"));

    G4cout << " ExcitE be4 Correct : " <<GetExcitationEnergy() << G4endl;
    G4cout << " Mom Transfered to nucleus : " << theMomentumTransfer << " " << theMomentumTransfer.mag() << G4endl;
    PrintKTVector(&theFinalState,std::string(" FinalState uncorrected"));
#endif

    //


    G4double ExcitationEnergy=GetExcitationEnergy();

#ifdef debug_BIC_Propagate_finals
    PrintKTVector(&theFinalState,std::string(" FinalState be4 corr"));
    G4cout << " Excitation Energy prefinal,  #collisions:, out, captured  "
            << ExcitationEnergy << " "
            << collisionCount << " "
            << theFinalState.size() << " "
            << theCapturedList.size()<<G4endl;
#endif

    if (ExcitationEnergy < 0 )
    {
        G4int maxtry=5, ntry=0;
        do {
            CorrectFinalPandE();
            ExcitationEnergy=GetExcitationEnergy();
        } while ( ++ntry < maxtry && ExcitationEnergy < 0 );       /* Loop checking, 31.08.2015, G.Folger */
    }
    _DebugEpConservation("corrected");

#ifdef debug_BIC_Propagate_finals
    PrintKTVector(&theFinalState,std::string(" FinalState corrected"));
    G4cout << " Excitation Energy final,  #collisions:, out, captured  "
            << ExcitationEnergy << " "
            << collisionCount << " "
            << theFinalState.size() << " "
            << theCapturedList.size()<<G4endl;
#endif


    if ( ExcitationEnergy < 0. )
    {
            #ifdef debug_G4BinaryCascade
              	  G4cerr << "G4BinaryCascade-Warning: negative excitation energy ";
              	  G4cerr <<ExcitationEnergy<<G4endl;
             	   PrintKTVector(&theFinalState,std::string("FinalState"));
             	   PrintKTVector(&theCapturedList,std::string("captured"));
             	  G4cout << "negative ExE:Final 4Momentum .mag: " << GetFinal4Momentum()
             	          << " "<< GetFinal4Momentum().mag()<< G4endl
             	          << "negative ExE:FinalNucleusMom  .mag: " << GetFinalNucleusMomentum()
             		  << " "<< GetFinalNucleusMomentum().mag()<< G4endl;
            #endif
            #ifdef debug_BIC_return
                    G4cout << "  negative Excitation E return empty products " << products << G4endl;
                    G4cout << "return 4, excit < 0 "<< G4endl;
            #endif

        ClearAndDestroy(products);
        return products;   // return empty products- FixMe
    }

    G4ReactionProductVector * precompoundProducts=DeExcite();


    G4DecayKineticTracks decay(&theFinalState);
    _DebugEpConservation("decayed");

    products= ProductsAddFinalState(products, theFinalState);

    products= ProductsAddPrecompound(products, precompoundProducts);

//    products=ProductsAddFakeGamma(products);


    thePrimaryEscape = true;

    #ifdef debug_BIC_return
    G4cout << "BIC: return @end, all ok "<< G4endl;
    //G4cout << "  return products " << products << G4endl;
    #endif

    return products;
}

//----------------------------------------------------------------------------
G4double G4BinaryCascade::GetExcitationEnergy()
//----------------------------------------------------------------------------
{

    // get A and Z for the residual nucleus
#if defined(debug_G4BinaryCascade) || defined(debug_BIC_GetExcitationEnergy)
    G4int finalA = theTargetList.size()+theCapturedList.size();
    G4int finalZ = GetTotalCharge(theTargetList)+GetTotalCharge(theCapturedList);
    if ( (currentA - finalA) != 0 || (currentZ - finalZ) != 0 )
    {
        G4cerr << "G4BIC:GetExcitationEnergy(): Nucleon counting error current/final{A,Z} "
                << "("<< currentA << "," << finalA << ") ("<< currentZ << "," << finalZ << ")" << G4endl;
    }

#endif

    G4double excitationE(0);
    G4double nucleusMass(0);
    if(currentZ>.5)
    {
        nucleusMass = GetIonMass(currentZ,currentA);
    }
    else if (currentZ==0 )     // Uzhi && currentA==1 )                     // Uzhi
    {                                                                       // Uzhi
        if(currentA == 1) {nucleusMass = G4Neutron::Neutron()->GetPDGMass();}// Uzhi
        else              {nucleusMass = GetFinalNucleusMomentum().mag()     // Uzhi
            - 3.*MeV*currentA;}                 // Uzhi
    }                                                                       // Uzhi
    else
    {
#ifdef debug_G4BinaryCascade
        G4cout << "G4BinaryCascade::GetExcitationEnergy(): Warning - invalid nucleus (A,Z)=("
                << currentA << "," << currentZ << ")" << G4endl;
#endif
        return 0;
    }

#ifdef debug_BIC_GetExcitationEnergy
    G4ping debug("debug_ExcitationEnergy");
    debug.push_back("====> current A, Z");
    debug.push_back(currentZ);
    debug.push_back(currentA);
    debug.push_back("====> final A, Z");
    debug.push_back(finalZ);
    debug.push_back(finalA);
    debug.push_back(nucleusMass);
    debug.push_back(GetFinalNucleusMomentum().mag());
    debug.dump();
    //  PrintKTVector(&theTargetList, std::string(" current target list info"));
    //PrintKTVector(&theCapturedList, std::string(" current captured list info"));
#endif

    excitationE = GetFinalNucleusMomentum().mag() - nucleusMass;

    //G4double exE2 = GetFinal4Momentum().mag() - nucleusMass;

    //G4cout << "old/new excitE " << excitationE << " / "<< exE2 << G4endl;

#ifdef debug_BIC_GetExcitationEnergy
    // ------ debug
    if ( excitationE < 0 )
    {
        G4cout << "negative ExE final Ion mass " <<nucleusMass<< G4endl;
        G4LorentzVector Nucl_mom=GetFinalNucleusMomentum();
        if(finalZ>.5) G4cout << " Final nuclmom/mass " << Nucl_mom << " " << Nucl_mom.mag()
		    		           << " (A,Z)=("<< finalA <<","<<finalZ <<")"
		    		           << " mass " << nucleusMass << " "
		    		           << " excitE " << excitationE << G4endl;


        G4int A = the3DNucleus->GetMassNumber();
        G4int Z = the3DNucleus->GetCharge();
        G4double initialExc(0);
        if(Z>.5)
        {
            initialExc = theInitial4Mom.mag()- GetIonMass(Z, A);
            G4cout << "GetExcitationEnergy: Initial nucleus A Z " << A << " " << Z << " " << initialExc << G4endl;
        }
    }

#endif

    return excitationE;
}


//----------------------------------------------------------------------------
//
//       P R I V A T E   M E M B E R   F U N C T I O N S
//
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
void G4BinaryCascade::BuildTargetList()
//----------------------------------------------------------------------------
{

    if(!the3DNucleus->StartLoop())
    {
        //    G4cerr << "G4BinaryCascade::BuildTargetList(): StartLoop() error!"
        //	   << G4endl;
        return;
    }

    ClearAndDestroy(&theTargetList);  // clear theTargetList before rebuilding

    G4Nucleon * nucleon;
    const G4ParticleDefinition * definition;
    G4ThreeVector pos;
    G4LorentzVector mom;
    // if there are nucleon hit by higher energy models, then SUM(momenta) != 0
    initialZ=the3DNucleus->GetCharge();
    initialA=the3DNucleus->GetMassNumber();
    initial_nuclear_mass=GetIonMass(initialZ,initialA);
    theInitial4Mom = G4LorentzVector(0,0,0,initial_nuclear_mass);
    currentA=0;
    currentZ=0;
    while((nucleon = the3DNucleus->GetNextNucleon()) != NULL)       /* Loop checking, 31.08.2015, G.Folger */
    {
        // check if nucleon is hit by higher energy model.
        if ( ! nucleon->AreYouHit() )
        {
            definition = nucleon->GetDefinition();
            pos = nucleon->GetPosition();
            mom = nucleon->GetMomentum();
            //    G4cout << "Nucleus " << pos.mag()/fermi << " " << mom.e() << G4endl;
            //theInitial4Mom += mom;
            //        the potential inside the nucleus is taken into account, and nucleons are on mass shell.
            mom.setE( std::sqrt( mom.vect().mag2() + sqr(definition->GetPDGMass()) ) );
            G4KineticTrack * kt = new G4KineticTrack(definition, 0., pos, mom);
            kt->SetState(G4KineticTrack::inside);
            kt->SetNucleon(nucleon);
            theTargetList.push_back(kt);
            ++currentA;
            if (definition->GetPDGCharge() > .5 ) ++currentZ;
        }

#ifdef debug_BIC_BuildTargetList
        else { G4cout << "nucleon is hit" << nucleon << G4endl;}
#endif
    }
    massInNucleus = 0;
    if(currentZ>.5)
    {
        massInNucleus = GetIonMass(currentZ,currentA);
    } else if (currentZ==0 && currentA>=1 )
    {
        massInNucleus = currentA * G4Neutron::Neutron()->GetPDGMass();
    } else
    {
        G4cerr << "G4BinaryCascade::BuildTargetList(): Fatal Error - invalid nucleus (A,Z)=("
                << currentA << "," << currentZ << ")" << G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde::BuildTargetList()");
    }
    currentInitialEnergy=	theInitial4Mom.e() + theProjectile4Momentum.e();

#ifdef debug_BIC_BuildTargetList
    G4cout << "G4BinaryCascade::BuildTargetList():  nucleus (A,Z)=("
            << currentA << "," << currentZ << ") mass: " << massInNucleus <<
            ", theInitial4Mom " << theInitial4Mom <<
            ", currentInitialEnergy " << currentInitialEnergy << G4endl;
#endif

}

//----------------------------------------------------------------------------
G4bool  G4BinaryCascade::BuildLateParticleCollisions(G4KineticTrackVector * secondaries)
//----------------------------------------------------------------------------
{
   G4bool success(false);
   std::vector<G4KineticTrack *>::iterator iter;

   lateA=lateZ=0;
   projectileA=projectileZ=0;

   G4double StartingTime=DBL_MAX;        // Search for minimal formation time
   for(iter = secondaries->begin(); iter != secondaries->end(); ++iter)
   {
      if((*iter)->GetFormationTime() < StartingTime)
         StartingTime = (*iter)->GetFormationTime();
   }

   //PrintKTVector(secondaries, "initial late particles ");
   G4LorentzVector lateParticles4Momentum(0,0,0,0);
   for(iter = secondaries->begin(); iter != secondaries->end(); ++iter)
   {
      //  G4cout << " Formation time : " << (*iter)->GetDefinition()->GetParticleName() << " "
      //   << (*iter)->GetFormationTime() << G4endl;
      G4double FormTime = (*iter)->GetFormationTime() - StartingTime;
      (*iter)->SetFormationTime(FormTime);
      if( (*iter)->GetState() == G4KineticTrack::undefined  )   // particles from high energy generator
      {
         FindLateParticleCollision(*iter);
         lateParticles4Momentum += (*iter)->GetTrackingMomentum();
         lateA += (*iter)->GetDefinition()->GetBaryonNumber();
         lateZ += G4lrint((*iter)->GetDefinition()->GetPDGCharge()/eplus);
         //PrintKTVector(*iter, "late particle ");
      } else
      {
         theSecondaryList.push_back(*iter);
         //PrintKTVector(*iter, "incoming particle ");
         theProjectile4Momentum += (*iter)->GetTrackingMomentum();
         projectileA += (*iter)->GetDefinition()->GetBaryonNumber();
         projectileZ += G4lrint((*iter)->GetDefinition()->GetPDGCharge()/eplus);
#ifdef debug_BIC_Propagate
         G4cout << " Adding initial secondary " << *iter
               << " time" << (*iter)->GetFormationTime()
               << ", state " << (*iter)->GetState() << G4endl;
#endif
      }
   }
   //theCollisionMgr->Print();
   const G4HadProjectile * primary = GetPrimaryProjectile();  // check for primary from TheoHE model

   if (primary){
      G4LorentzVector mom=primary->Get4Momentum();
      theProjectile4Momentum += mom;
      projectileA = primary->GetDefinition()->GetBaryonNumber();
      projectileZ = G4lrint(primary->GetDefinition()->GetPDGCharge()/eplus);
      // now check if "excitation" energy left by TheoHE model
      G4double excitation= theProjectile4Momentum.e() + initial_nuclear_mass - lateParticles4Momentum.e() - massInNucleus;
#ifdef debug_BIC_GetExcitationEnergy
      G4cout << "BIC: Proj.e, nucl initial, nucl final, lateParticles"
            << theProjectile4Momentum << ",  "
            << initial_nuclear_mass<< ",  " << massInNucleus << ",  "
            << lateParticles4Momentum << G4endl;
      G4cout << "BIC: Proj.e / initial excitation: " << theProjectile4Momentum.e() << " / " << excitation << G4endl;
#endif
      success = excitation > 0;
#ifdef debug_G4BinaryCascade
      if ( ! success ) {
         G4cout << "G4BinaryCascade::BuildLateParticleCollisions(): Proj.e / initial excitation: " << theProjectile4Momentum.e() << " / " << excitation << G4endl;
         //PrintKTVector(secondaries);
      }
#endif
   } else {
      // no primary from HE model -> cascade
      success=true;
   }

   if (success) {
      secondaries->clear(); // Don't leave "G4KineticTrack *"s in two vectors
      delete secondaries;
   }
   return success;
}

//----------------------------------------------------------------------------
G4ReactionProductVector *  G4BinaryCascade::DeExcite()
//----------------------------------------------------------------------------
{
   // find a fragment and call the precompound model.
   G4Fragment * fragment = 0;
   G4ReactionProductVector * precompoundProducts = 0;

   G4LorentzVector pFragment(0);
   // G4cout << " final4mon " << GetFinal4Momentum() /MeV << G4endl;

   //   if ( ExcitationEnergy >= 0 )                                         // closed by Uzhi
   //   {                                                                    // closed by Uzhi
   fragment = FindFragments();


   if(fragment)                                                            // Uzhi
   {                                                                       // Uzhi
      if(fragment->GetA_asInt() >1)                                          // Uzhi
      {
         pFragment=fragment->GetMomentum();
         // G4cout << " going to preco with fragment 4 mom " << pFragment << G4endl;
         if (theDeExcitation)                // pre-compound
         {
            precompoundProducts= theDeExcitation->DeExcite(*fragment);
         }
         else if (theExcitationHandler)    // de-excitation
         {
            precompoundProducts=theExcitationHandler->BreakItUp(*fragment);
         }

      } else
      {                                   // fragment->GetA_asInt() <= 1, so a single proton, as a fragment must have Z>0
         if (theTargetList.size() + theCapturedList.size() > 1 ) {
            throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde:: Invalid Fragment");
         }

         std::vector<G4KineticTrack *>::iterator i;
         if ( theTargetList.size() == 1 )  {i=theTargetList.begin();}
         if ( theCapturedList.size() == 1 ) {i=theCapturedList.begin();}                             // Uzhi
         G4ReactionProduct * aNew = new G4ReactionProduct((*i)->GetDefinition());
         aNew->SetTotalEnergy((*i)->GetDefinition()->GetPDGMass());
         aNew->SetCreatorModel(theBIC_ID);
         aNew->SetMomentum(G4ThreeVector(0));// see boost for preCompoundProducts below..
         precompoundProducts = new G4ReactionProductVector();
         precompoundProducts->push_back(aNew);
      }                            // End of fragment->GetA() < 1.5
      delete fragment;
      fragment=0;

   } else                            // End of if(fragment)
   {                                 // No fragment, can be neutrons only  // Uzhi

      precompoundProducts = DecayVoidNucleus();
   }
   #ifdef debug_BIC_DeexcitationProducts

       G4LorentzVector fragment_momentum=GetFinalNucleusMomentum();
       G4LorentzVector Preco_momentum;
       if ( precompoundProducts )
       {
          std::vector<G4ReactionProduct *>::iterator j;
          for(j = precompoundProducts->begin(); j != precompoundProducts->end(); ++j)
          {
             G4LorentzVector pProduct((*j)->GetMomentum(),(*j)->GetTotalEnergy());
             Preco_momentum += pProduct;
           }
       }
       G4cout << "finalNuclMom / sum preco products" << fragment_momentum << "  / " << Preco_momentum
    		   << " delta E "<< fragment_momentum.e() - Preco_momentum.e() <<  G4endl;

   #endif

   return precompoundProducts;
}

//----------------------------------------------------------------------------
G4ReactionProductVector *  G4BinaryCascade::DecayVoidNucleus()
//----------------------------------------------------------------------------
{
   G4ReactionProductVector * result=0;
   if ( (theTargetList.size()+theCapturedList.size()) > 0 )
   {
      result = new G4ReactionProductVector;
      std::vector<G4KineticTrack *>::iterator aNuc;
      G4LorentzVector aVec;
      std::vector<G4double> masses;
      G4double sumMass(0);

      if ( theTargetList.size() != 0)                                      // Uzhi
            {
         for ( aNuc=theTargetList.begin(); aNuc != theTargetList.end(); aNuc++)
         {
            G4double mass=(*aNuc)->GetDefinition()->GetPDGMass();
            masses.push_back(mass);
            sumMass += mass;
         }
            }                                                                    // Uzhi

      if ( theCapturedList.size() != 0)                                    // Uzhi
      {                                                                    // Uzhi
         for(aNuc = theCapturedList.begin();                               // Uzhi
               aNuc != theCapturedList.end(); aNuc++)                        // Uzhi
         {                                                                 // Uzhi
            G4double mass=(*aNuc)->GetDefinition()->GetPDGMass();          // Uzhi
            masses.push_back(mass);                                        // Uzhi
            sumMass += mass;                                               // Uzhi
         }
      }

      G4LorentzVector finalP=GetFinal4Momentum();
      G4FermiPhaseSpaceDecay decay;
      // G4cout << " some neutrons? " << masses.size() <<" " ;
      // G4cout<< theTargetList.size()<<" "<<finalP <<" " << finalP.mag()<<G4endl;

      G4double eCMS=finalP.mag();
      if ( eCMS < sumMass )                    // @@GF --- Cheat!!
      {
         eCMS=sumMass + 2*MeV*masses.size();
         finalP.setE(std::sqrt(finalP.vect().mag2() + sqr(eCMS)));
      }

      precompoundLorentzboost.set(finalP.boostVector());
      std::vector<G4LorentzVector*> * momenta=decay.Decay(eCMS,masses);
      std::vector<G4LorentzVector*>::iterator aMom=momenta->begin();


      if ( theTargetList.size() != 0)
      {
         for ( aNuc=theTargetList.begin();
               (aNuc != theTargetList.end()) && (aMom!=momenta->end());
               aNuc++, aMom++ )
         {
            G4ReactionProduct * aNew = new G4ReactionProduct((*aNuc)->GetDefinition());
            aNew->SetTotalEnergy((*aMom)->e());
            aNew->SetMomentum((*aMom)->vect());
            aNew->SetCreatorModel(theBIC_ID);
            result->push_back(aNew);

            delete *aMom;
         }
      }

      if ( theCapturedList.size() != 0)                                    // Uzhi
            {                                                                    // Uzhi
         for ( aNuc=theCapturedList.begin();                                // Uzhi
               (aNuc != theCapturedList.end()) && (aMom!=momenta->end());    // Uzhi
               aNuc++, aMom++ )                                             // Uzhi
         {                                                                  // Uzhi
            G4ReactionProduct * aNew = new G4ReactionProduct(               // Uzhi
                  (*aNuc)->GetDefinition());            // Uzhi
            aNew->SetTotalEnergy((*aMom)->e());                             // Uzhi
            aNew->SetMomentum((*aMom)->vect());                             // Uzhi
            aNew->SetCreatorModel(theBIC_ID);
            result->push_back(aNew);                           // Uzhi
            delete *aMom;                                                   // Uzhi
         }                                                                  // Uzhi
            }                                                                    // Uzhi

      delete momenta;
   }
   return result;
}                   // End if(!fragment)

//----------------------------------------------------------------------------
G4ReactionProductVector * G4BinaryCascade::ProductsAddFinalState(G4ReactionProductVector * products, G4KineticTrackVector & fs)
//----------------------------------------------------------------------------
{
// fill in products the outgoing particles
    size_t i(0);
#ifdef debug_BIC_Propagate_finals
    G4LorentzVector mom_fs;
#endif
    for(i = 0; i< fs.size(); i++)
    {
        G4KineticTrack * kt = fs[i];
        G4ReactionProduct * aNew = new G4ReactionProduct(kt->GetDefinition());
        aNew->SetMomentum(kt->Get4Momentum().vect());
        aNew->SetTotalEnergy(kt->Get4Momentum().e());
        aNew->SetNewlyAdded(kt->IsParticipant());
        aNew->SetCreatorModel(theBIC_ID);
        products->push_back(aNew);

#ifdef debug_BIC_Propagate_finals
        mom_fs += kt->Get4Momentum();
        G4cout <<kt->GetDefinition()->GetParticleName();
        G4cout << " Particle Ekin " << aNew->GetKineticEnergy();
        G4cout << ", is " << (kt->GetDefinition()->GetPDGStable() ? "stable" :
                (kt->GetDefinition()->IsShortLived() ? "short lived " : "non stable"))  ;
        G4cout << G4endl;
#endif

    }
#ifdef debug_BIC_Propagate_finals
    G4cout << " Final state momentum " << mom_fs << G4endl;
#endif

    return products;
}
//----------------------------------------------------------------------------
G4ReactionProductVector * G4BinaryCascade::ProductsAddPrecompound(G4ReactionProductVector * products, G4ReactionProductVector * precompoundProducts)
//----------------------------------------------------------------------------
{
   G4LorentzVector pSumPreco(0), pPreco(0);

   if ( precompoundProducts )
   {
      std::vector<G4ReactionProduct *>::iterator j;
      for(j = precompoundProducts->begin(); j != precompoundProducts->end(); ++j)
      {
         // boost back to system of moving nucleus
         G4LorentzVector pProduct((*j)->GetMomentum(),(*j)->GetTotalEnergy());
         pPreco+= pProduct;
#ifdef debug_BIC_Propagate_finals
         G4cout << "BIC: pProduct be4 boost " <<pProduct << G4endl;
#endif
         pProduct *= precompoundLorentzboost;
#ifdef debug_BIC_Propagate_finals
         G4cout << "BIC: pProduct aft boost " <<pProduct << G4endl;
#endif
         pSumPreco += pProduct;
         (*j)->SetTotalEnergy(pProduct.e());
         (*j)->SetMomentum(pProduct.vect());
         (*j)->SetNewlyAdded(true);
         products->push_back(*j);
      }
      // G4cout << " unboosted preco result mom " << pPreco / MeV << "  ..- fragmentMom " << (pPreco - pFragment)/MeV<< G4endl;
      // G4cout << " preco result mom " << pSumPreco / MeV << "  ..-file4Mom " << (pSumPreco - GetFinal4Momentum())/MeV<< G4endl;
      precompoundProducts->clear();
      delete precompoundProducts;
   }
   return products;
}
//----------------------------------------------------------------------------
void  G4BinaryCascade::FindCollisions(G4KineticTrackVector * secondaries)
//----------------------------------------------------------------------------
{
    for(std::vector<G4KineticTrack *>::iterator i = secondaries->begin();
            i != secondaries->end(); ++i)
    {
        for(std::vector<G4BCAction *>::iterator j = theImR.begin();
                j!=theImR.end(); j++)
        {
            //      G4cout << "G4BinaryCascade::FindCollisions: at action " << *j << G4endl;
            const std::vector<G4CollisionInitialState *> & aCandList
            = (*j)->GetCollisions(*i, theTargetList, theCurrentTime);
            for(size_t count=0; count<aCandList.size(); count++)
            {
                theCollisionMgr->AddCollision(aCandList[count]);
                //4cout << "====================== New Collision ================="<<G4endl;
                //theCollisionMgr->Print();
            }
        }
    }
}


//----------------------------------------------------------------------------
void  G4BinaryCascade::FindDecayCollision(G4KineticTrack * secondary)
//----------------------------------------------------------------------------
{
    const std::vector<G4CollisionInitialState *> & aCandList
    = theDecay->GetCollisions(secondary, theTargetList, theCurrentTime);
    for(size_t count=0; count<aCandList.size(); count++)
    {
        theCollisionMgr->AddCollision(aCandList[count]);
    }
}

//----------------------------------------------------------------------------
void  G4BinaryCascade::FindLateParticleCollision(G4KineticTrack * secondary)
//----------------------------------------------------------------------------
{

    G4double tin=0., tout=0.;
    if (((G4RKPropagation*)thePropagator)->GetSphereIntersectionTimes(secondary,tin,tout))
    {
        if ( tin > 0 )
        {
            secondary->SetState(G4KineticTrack::outside);
        } else if ( tout > 0 )
        {
            secondary->SetState(G4KineticTrack::inside);
        } else {
            //G4cout << "G4BC set miss , tin, tout " << tin << " , " << tout <<G4endl;
            secondary->SetState(G4KineticTrack::miss_nucleus);
        }
    } else {
        secondary->SetState(G4KineticTrack::miss_nucleus);
        //G4cout << "G4BC set miss ,no intersect tin, tout " << tin << " , " << tout <<G4endl;
    }


#ifdef debug_BIC_FindCollision
    G4cout << "FindLateP Particle, 4-mom, times newState "
            << secondary->GetDefinition()->GetParticleName() << " "
            << secondary->Get4Momentum()
            << " times " <<  tin << " " << tout << " "
            << secondary->GetState() << G4endl;
#endif

    const std::vector<G4CollisionInitialState *> & aCandList
    = theLateParticle->GetCollisions(secondary, theTargetList, theCurrentTime);
    for(size_t count=0; count<aCandList.size(); count++)
    {
#ifdef debug_BIC_FindCollision
        G4cout << " Adding a late Col : " << aCandList[count] << G4endl;
#endif
        theCollisionMgr->AddCollision(aCandList[count]);
    }
}


//----------------------------------------------------------------------------
G4bool G4BinaryCascade::ApplyCollision(G4CollisionInitialState * collision)
//----------------------------------------------------------------------------
{
    G4KineticTrack * primary = collision->GetPrimary();

#ifdef debug_BIC_ApplyCollision
    G4cerr << "G4BinaryCascade::ApplyCollision start"<<G4endl;
    theCollisionMgr->Print();
    G4cout << "ApplyCollisions : projte 4mom " << primary->GetTrackingMomentum()<< G4endl;
#endif

    G4KineticTrackVector target_collection=collision->GetTargetCollection();
    G4bool haveTarget=target_collection.size()>0;
    if( haveTarget && (primary->GetState() != G4KineticTrack::inside) )
    {
#ifdef debug_G4BinaryCascade
        G4cout << "G4BinaryCasacde::ApplyCollision(): StateError " << primary << G4endl;
        PrintKTVector(primary,std::string("primay- ..."));
        PrintKTVector(&target_collection,std::string("... targets"));
        collision->Print();
        G4cout << G4endl;
        theCollisionMgr->Print();
        //*GF*     throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde::ApplyCollision()");
#endif
        return false;
//    } else {
//       G4cout << "G4BinaryCasacde::ApplyCollision(): decay " << G4endl;
//       PrintKTVector(primary,std::string("primay- ..."));
//       G4double tin=0., tout=0.;
//       if (((G4RKPropagation*)thePropagator)->GetSphereIntersectionTimes(primary,tin,tout))
//       {
//           G4cout << "tin tout: " << tin << " " << tout << G4endl;
//       }

    }

    G4LorentzVector mom4Primary=primary->Get4Momentum();

    G4int initialBaryon(0);
    G4int initialCharge(0);
    if (primary->GetState() == G4KineticTrack::inside)
    {
        initialBaryon = primary->GetDefinition()->GetBaryonNumber();
        initialCharge = G4lrint(primary->GetDefinition()->GetPDGCharge()/eplus);
    }

    // for primary resonances, subtract neutron ( = proton) field ( ie. add std::abs(field))
    G4double initial_Efermi=CorrectShortlivedPrimaryForFermi(primary,target_collection);
    //****************************************


    G4KineticTrackVector * products = collision->GetFinalState();

#ifdef debug_BIC_ApplyCollision
    DebugApplyCollisionFail(collision, products);
#endif

    // reset primary to initial state, in case there is a veto...
    primary->Set4Momentum(mom4Primary);

    G4bool lateParticleCollision= (!haveTarget) && products && products->size() == 1;
    G4bool decayCollision= (!haveTarget) && products && products->size() > 1;
    G4bool Success(true);


#ifdef debug_G4BinaryCascade
    G4int lateBaryon(0), lateCharge(0);
#endif

    if ( lateParticleCollision )
    {  // for late particles, reset charges
        //G4cout << "lateP, initial B C state " << initialBaryon << " "
        //        << initialCharge<< " " << primary->GetState() << " "<< primary->GetDefinition()->GetParticleName()<< G4endl;
#ifdef debug_G4BinaryCascade
        lateBaryon = initialBaryon;
        lateCharge = initialCharge;
#endif
        initialBaryon=initialCharge=0;
        lateA -= primary->GetDefinition()->GetBaryonNumber();
        lateZ -= G4lrint(primary->GetDefinition()->GetPDGCharge()/eplus);
    }

    initialBaryon += collision->GetTargetBaryonNumber();
    initialCharge += G4lrint(collision->GetTargetCharge());
    if (!lateParticleCollision)
    {
       if( !products || products->size()==0 || !CheckPauliPrinciple(products) )
       {
#ifdef debug_BIC_ApplyCollision
          if (products) G4cout << " ======Failed Pauli =====" << G4endl;
          G4cerr << "G4BinaryCascade::ApplyCollision blocked"<<G4endl;
#endif
          Success=false;
       }



       if (Success && primary->GetState() == G4KineticTrack::inside ) {   // if the primary was outside, nothing to correct
          if (! CorrectShortlivedFinalsForFermi(products, initial_Efermi)){
             Success=false;
          }
       }
    }

#ifdef debug_BIC_ApplyCollision
    DebugApplyCollision(collision, products);
#endif

    if ( ! Success ){
        if (products) ClearAndDestroy(products);
        if ( decayCollision ) FindDecayCollision(primary);  // for decay, sample new decay
        delete products;
        products=0;
        return false;
    }

    G4int finalBaryon(0);
    G4int finalCharge(0);
    G4KineticTrackVector toFinalState;
    for(std::vector<G4KineticTrack *>::iterator i =products->begin(); i != products->end(); i++)
    {
        if ( ! lateParticleCollision )
        {
            (*i)->SetState(primary->GetState());  // decay may be anywhere!
            if ( (*i)->GetState() == G4KineticTrack::inside ){
                finalBaryon+=(*i)->GetDefinition()->GetBaryonNumber();
                finalCharge+=G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus);
            } else {
               G4double tin=0., tout=0.;
               if (((G4RKPropagation*)thePropagator)->GetSphereIntersectionTimes((*i),tin,tout) &&
                     tin < 0 && tout > 0 )
               {
                  PrintKTVector((*i),"particle inside marked not-inside");
                   G4cout << "tin tout: " << tin << " " << tout << G4endl;
               }
            }
        } else {
            G4double tin=0., tout=0.;
            if (((G4RKPropagation*)thePropagator)->GetSphereIntersectionTimes((*i),tin,tout))
            {
                //G4cout << "tin tout: " << tin << " " << tout << G4endl;
                if ( tin > 0 )
                {
                    (*i)->SetState(G4KineticTrack::outside);
                }
                else if ( tout > 0 )
                {
                    (*i)->SetState(G4KineticTrack::inside);
                    finalBaryon+=(*i)->GetDefinition()->GetBaryonNumber();
                    finalCharge+=G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus);
                }
                else
                {
                    (*i)->SetState(G4KineticTrack::gone_out);
                    toFinalState.push_back((*i));
                }
            } else
            {
                (*i)->SetState(G4KineticTrack::miss_nucleus);
                //G4cout << " G4BC - miss -late Part- no intersection found " << G4endl;
                toFinalState.push_back((*i));
            }

        }
    }
    if(!toFinalState.empty())
    {
        theFinalState.insert(theFinalState.end(),
                toFinalState.begin(),toFinalState.end());
        std::vector<G4KineticTrack *>::iterator iter1, iter2;
        for(iter1 = toFinalState.begin(); iter1 != toFinalState.end();
                ++iter1)
        {
            iter2 = std::find(products->begin(), products->end(),
                    *iter1);
            if ( iter2 != products->end() ) products->erase(iter2);
        }
        theCollisionMgr->RemoveTracksCollisions(&toFinalState);
    }

    //G4cout << " currentA, Z be4: " << currentA << " " << currentZ << G4endl;
    currentA += finalBaryon-initialBaryon;
    currentZ += finalCharge-initialCharge;
    //G4cout << " ApplyCollision currentA, Z aft: " << currentA << " " << currentZ << G4endl;

    G4KineticTrackVector oldSecondaries;
    oldSecondaries.push_back(primary);
    primary->Hit();

#ifdef debug_G4BinaryCascade
    if ( (finalBaryon-initialBaryon-lateBaryon) != 0 || (finalCharge-initialCharge-lateCharge) != 0 )
    {
        G4cout << "G4BinaryCascade: Error in Balancing: " << G4endl;
        G4cout << "initial/final baryon number, initial/final Charge "
                << initialBaryon <<" "<< finalBaryon <<" "
                << initialCharge <<" "<< finalCharge <<" "
                << " in Collision type: "<< typeid(*collision->GetGenerator()).name()
                << ", with number of products: "<< products->size() <<G4endl;
        G4cout << G4endl<<"Initial condition are these:"<<G4endl;
        G4cout << "proj: "<<collision->GetPrimary()->GetDefinition()->GetParticleName()<<G4endl;
        for(size_t it=0; it<collision->GetTargetCollection().size(); it++)
        {
            G4cout << "targ: "
                    <<collision->GetTargetCollection()[it]->GetDefinition()->GetParticleName()<<G4endl;
        }
        PrintKTVector(&collision->GetTargetCollection(),std::string(" Target particles"));
        G4cout << G4endl<<G4endl;
    }
#endif

    G4KineticTrackVector oldTarget = collision->GetTargetCollection();
    for(size_t ii=0; ii< oldTarget.size(); ii++)
    {
        oldTarget[ii]->Hit();
    }

    UpdateTracksAndCollisions(&oldSecondaries, &oldTarget, products);
    std::for_each(oldSecondaries.begin(), oldSecondaries.end(), Delete<G4KineticTrack>());
    std::for_each(oldTarget.begin(), oldTarget.end(), Delete<G4KineticTrack>());

    delete products;
    return true;
}


//----------------------------------------------------------------------------
G4bool G4BinaryCascade::Absorb()
//----------------------------------------------------------------------------
{
    // Do it in two step: cannot change theSecondaryList inside the first loop.
    G4Absorber absorber(theCutOnPAbsorb);

    // Build the vector of G4KineticTracks that must be absorbed
    G4KineticTrackVector absorbList;
    std::vector<G4KineticTrack *>::iterator iter;
    //  PrintKTVector(&theSecondaryList, " testing for Absorb" );
    for(iter = theSecondaryList.begin();
            iter != theSecondaryList.end(); ++iter)
    {
        G4KineticTrack * kt = *iter;
        if(kt->GetState() == G4KineticTrack::inside)// absorption happens only inside the nucleus
        {
            if(absorber.WillBeAbsorbed(*kt))
            {
                absorbList.push_back(kt);
            }
        }
    }

    if(absorbList.empty())
        return false;

    G4KineticTrackVector toDelete;
    for(iter = absorbList.begin(); iter != absorbList.end(); ++iter)
    {
        G4KineticTrack * kt = *iter;
        if(!absorber.FindAbsorbers(*kt, theTargetList))
            throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::Absorb(): Cannot absorb a particle.");

        if(!absorber.FindProducts(*kt))
            throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::Absorb(): Cannot absorb a particle.");

        G4KineticTrackVector * products = absorber.GetProducts();
        G4int maxLoopCount = 1000;
        while(!CheckPauliPrinciple(products) && --maxLoopCount>0)   /* Loop checking, 31.08.2015, G.Folger */
        {
            ClearAndDestroy(products);
            if(!absorber.FindProducts(*kt))
                throw G4HadronicException(__FILE__, __LINE__,
                        "G4BinaryCascade::Absorb(): Cannot absorb a particle.");
        }
		  if ( --maxLoopCount < 0 ) throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::Absorb(): Cannot absorb a particle."); 
        // ------ debug
        //    G4cerr << "Absorb CheckPauliPrinciple count= " <<  maxLoopCount << G4endl;
        // ------ end debug
        G4KineticTrackVector toRemove;  // build a vector for UpdateTrack...
        toRemove.push_back(kt);
        toDelete.push_back(kt);  // delete the track later
        G4KineticTrackVector * absorbers = absorber.GetAbsorbers();
        UpdateTracksAndCollisions(&toRemove, absorbers, products);
        ClearAndDestroy(absorbers);
    }
    ClearAndDestroy(&toDelete);
    return true;
}



// Capture all p and n with Energy < theCutOnP
//----------------------------------------------------------------------------
G4bool G4BinaryCascade::Capture(G4bool verbose)
//----------------------------------------------------------------------------
{
    G4KineticTrackVector captured;
    G4bool capture = false;
    std::vector<G4KineticTrack *>::iterator i;

    G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;

    G4double capturedEnergy = 0;
    G4int particlesAboveCut=0;
    G4int particlesBelowCut=0;
    if ( verbose ) G4cout << " Capture: secondaries " << theSecondaryList.size() << G4endl;
    for(i = theSecondaryList.begin(); i != theSecondaryList.end(); ++i)
    {
        G4KineticTrack * kt = *i;
        if (verbose) G4cout << "Capture position, radius, state " <<kt->GetPosition().mag()<<" "<<theOuterRadius<<" "<<kt->GetState()<<G4endl;
        if(kt->GetState() == G4KineticTrack::inside) // capture happens only inside the nucleus
        {
            if((kt->GetDefinition() == G4Proton::Proton()) ||
                    (kt->GetDefinition() == G4Neutron::Neutron()))
            {
                //GF cut on kinetic energy    if(kt->Get4Momentum().vect().mag() >= theCutOnP)
                G4double field=RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition())
	            		           -RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding());
                G4double energy= kt->Get4Momentum().e() - kt->GetActualMass() + field;
                if (verbose ) G4cout << "Capture: .e(), mass, field, energy" << kt->Get4Momentum().e() <<" "<<kt->GetActualMass()<<" "<<field<<" "<<energy<< G4endl;
                //	 if( energy < theCutOnP )
                //	 {
                capturedEnergy+=energy;
                ++particlesBelowCut;
                //	 } else
                //	 {
                //	    ++particlesAboveCut;
                //	 }
            }
        }
    }
    if (verbose) G4cout << "Capture particlesAboveCut,particlesBelowCut, capturedEnergy,capturedEnergy/particlesBelowCut <? 0.2*theCutOnP "
            << particlesAboveCut << " " << particlesBelowCut << " " << capturedEnergy
            << " "  << G4endl;
//    << " " << (particlesBelowCut>0) ? (capturedEnergy/particlesBelowCut) : (capturedEnergy) << " " << 0.2*theCutOnP << G4endl;
    //  if(particlesAboveCut==0 && particlesBelowCut>0 && capturedEnergy/particlesBelowCut<0.2*theCutOnP)
    if(particlesBelowCut>0 && capturedEnergy/particlesBelowCut<0.2*theCutOnP)
    {
        capture=true;
        for(i = theSecondaryList.begin(); i != theSecondaryList.end(); ++i)
        {
            G4KineticTrack * kt = *i;
            if(kt->GetState() == G4KineticTrack::inside) // capture happens only inside the nucleus
            {
                if((kt->GetDefinition() == G4Proton::Proton()) ||
                        (kt->GetDefinition() == G4Neutron::Neutron()))
                {
                    captured.push_back(kt);
                    kt->Hit();				//
                    theCapturedList.push_back(kt);
                }
            }
        }
        UpdateTracksAndCollisions(&captured, NULL, NULL);
    }

    return capture;
}


//----------------------------------------------------------------------------
G4bool G4BinaryCascade::CheckPauliPrinciple(G4KineticTrackVector * products)
//----------------------------------------------------------------------------
{
    G4int A = the3DNucleus->GetMassNumber();
    G4int Z = the3DNucleus->GetCharge();

    G4FermiMomentum fermiMom;
    fermiMom.Init(A, Z);

    const G4VNuclearDensity * density=the3DNucleus->GetNuclearDensity();

    G4KineticTrackVector::iterator i;
    const G4ParticleDefinition * definition;

    // ------ debug
    G4bool myflag = true;
    // ------ end debug
    //  G4ThreeVector xpos(0);
    for(i = products->begin(); i != products->end(); ++i)
    {
        definition = (*i)->GetDefinition();
        if((definition == G4Proton::Proton()) ||
                (definition == G4Neutron::Neutron()))
        {
            G4ThreeVector pos = (*i)->GetPosition();
            G4double d = density->GetDensity(pos);
            // energy correspondiing to fermi momentum
            G4double eFermi = std::sqrt( sqr(fermiMom.GetFermiMomentum(d)) + (*i)->Get4Momentum().mag2() );
            if( definition == G4Proton::Proton() )
            {
                eFermi -= the3DNucleus->CoulombBarrier();
            }
            G4LorentzVector mom = (*i)->Get4Momentum();
            // ------ debug
            /*
             *        G4cout << "p:[" << (1/MeV)*mom.x() << " " << (1/MeV)*mom.y() << " "
             *            << (1/MeV)*mom.z() << "] |p3|:"
             *            << (1/MeV)*mom.vect().mag() << " E:" << (1/MeV)*mom.t() << " ] m: "
             *            << (1/MeV)*mom.mag() << "  pos["
             *            << (1/fermi)*pos.x() << " "<< (1/fermi)*pos.y() << " "
             *            << (1/fermi)*pos.z() << "] |Dpos|: "
             *            << (1/fermi)*(pos-xpos).mag() << " Pfermi:"
             *            << (1/MeV)*p << G4endl;
             *         xpos=pos;
             */
            // ------ end debug
            if(mom.e() < eFermi )
            {
                // ------ debug
                myflag = false;
                // ------ end debug
                //      return false;
            }
        }
    }
#ifdef debug_BIC_CheckPauli
    if ( myflag  )
    {
        for(i = products->begin(); i != products->end(); ++i)
        {
            definition = (*i)->GetDefinition();
            if((definition == G4Proton::Proton()) ||
                    (definition == G4Neutron::Neutron()))
            {
                G4ThreeVector pos = (*i)->GetPosition();
                G4double d = density->GetDensity(pos);
                G4double pFermi = fermiMom.GetFermiMomentum(d);
                G4LorentzVector mom = (*i)->Get4Momentum();
                G4double field =((G4RKPropagation*)thePropagator)->GetField(definition->GetPDGEncoding(),pos);
                if ( mom.e()-mom.mag()+field > 160*MeV )
                {
                    G4cout << "momentum problem pFermi=" <<  pFermi
                            << " mom, mom.m " << mom << " " << mom.mag()
                            << " field " << field << G4endl;
                }
            }
        }
    }
#endif

    return myflag;
}

//----------------------------------------------------------------------------
void G4BinaryCascade::StepParticlesOut()
//----------------------------------------------------------------------------
{
    G4int counter=0;
    G4int countreset=0;
    //G4cout << " nucl. Radius " << radius << G4endl;
    // G4cerr <<"pre-while- theSecondaryList "<<G4endl;
    while( theSecondaryList.size() > 0 )               /* Loop checking, 31.08.2015, G.Folger */
	                                                    // if countreset reaches limit, there is a break from while, see below.
    {
        G4int nsec=0;
        G4double minTimeStep = 1.e-12*ns;   // about 30*fermi/(0.1*c_light);1.e-12*ns
                                            // i.e. a big step
        std::vector<G4KineticTrack *>::iterator i;
        for(i = theSecondaryList.begin(); i != theSecondaryList.end(); ++i)
        {
            G4KineticTrack * kt = *i;
            if( kt->GetState() == G4KineticTrack::inside )
            {
                nsec++;
                G4double tStep(0), tdummy(0);
                G4bool intersect =
                        ((G4RKPropagation*)thePropagator)->GetSphereIntersectionTimes(kt,tdummy,tStep);
#ifdef debug_BIC_StepParticlesOut
                G4cout << " minTimeStep, tStep Particle " <<minTimeStep << " " <<tStep
                        << " " <<kt->GetDefinition()->GetParticleName()
                        << " 4mom " << kt->GetTrackingMomentum()<<G4endl;
                if ( ! intersect );
                {
                    PrintKTVector(&theSecondaryList, std::string(" state ERROR....."));
                    throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::StepParticlesOut() particle not in nucleus");
                }
#endif
                if(intersect && tStep<minTimeStep && tStep> 0 )
                {
                    minTimeStep = tStep;
                }
            } else if ( kt->GetState() != G4KineticTrack::outside ){
                PrintKTVector(&theSecondaryList, std::string(" state ERROR....."));
                throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::StepParticlesOut() particle not in nucleus");
            }
        }
        minTimeStep *= 1.2;
        //    G4cerr << "CaptureCount = "<<counter<<" "<<nsec<<" "<<minTimeStep<<" "<<1*ns<<G4endl;
        G4double timeToCollision=DBL_MAX;
        G4CollisionInitialState * nextCollision=0;
        if(theCollisionMgr->Entries() > 0)
        {
            nextCollision = theCollisionMgr->GetNextCollision();
            timeToCollision = nextCollision->GetCollisionTime()-theCurrentTime;
            //  G4cout << " NextCollision  * , Time= " << nextCollision << " " <<timeToCollision<< G4endl;
        }
        if ( timeToCollision > minTimeStep )
        {
            DoTimeStep(minTimeStep);
            ++counter;
        } else
        {
            if (!DoTimeStep(timeToCollision) )
            {
                // Check if nextCollision is still valid, ie. partcile did not leave nucleus
                if (theCollisionMgr->GetNextCollision() != nextCollision )
                {
                    nextCollision = 0;
                }
            }
            // G4cerr <<"post- DoTimeStep 3"<<G4endl;

            if(nextCollision)
            {
                if  ( ApplyCollision(nextCollision))
                {
                    // G4cout << "ApplyCollision sucess " << G4endl;
                } else
                {
                    theCollisionMgr->RemoveCollision(nextCollision);
                }
            }
        }

        if(countreset>100)
        {
#ifdef debug_G4BinaryCascade
            G4cerr << "G4BinaryCascade.cc: Warning - aborting looping particle(s)" << G4endl;
            PrintKTVector(&theSecondaryList," looping particles added to theFinalState");
#endif

            //  add left secondaries to FinalSate
            std::vector<G4KineticTrack *>::iterator iter;
            for ( iter =theSecondaryList.begin(); iter != theSecondaryList.end(); ++iter)
            {
                theFinalState.push_back(*iter);
            }
            theSecondaryList.clear();

            break;
        }

        if(Absorb())
        {
            //       haveProducts = true;
            // G4cout << "Absorb sucess " << G4endl;
        }

        if(Capture(false))
        {
            //       haveProducts = true;
#ifdef debug_BIC_StepParticlesOut
            G4cout << "Capture sucess " << G4endl;
#endif
        }
        if ( counter > 100 && theCollisionMgr->Entries() == 0)   // no collision, and stepping for some time....
        {
#ifdef debug_BIC_StepParticlesOut
            PrintKTVector(&theSecondaryList,std::string("stepping 100 steps"));
#endif
            FindCollisions(&theSecondaryList);
            counter=0;
            ++countreset;
        }
        //G4cout << "currentZ @ end loop " << currentZ << G4endl;
        if ( ! currentZ ){
            // nucleus completely destroyed, fill in ReactionProductVector
           // products = FillVoidNucleusProducts(products);
    #ifdef debug_BIC_return
            G4cout << "return @ Z=0 after collision loop "<< G4endl;
            PrintKTVector(&theSecondaryList,std::string(" theSecondaryList"));
            G4cout << "theTargetList size: " << theTargetList.size() << G4endl;
            PrintKTVector(&theTargetList,std::string(" theTargetList"));
            PrintKTVector(&theCapturedList,std::string(" theCapturedList"));

            G4cout << " ExcitE be4 Correct : " <<GetExcitationEnergy() << G4endl;
            G4cout << " Mom Transfered to nucleus : " << theMomentumTransfer << " " << theMomentumTransfer.mag() << G4endl;
            PrintKTVector(&theFinalState,std::string(" FinalState uncorrected"));
          //  G4cout << "returned products: " << products->size() << G4endl;
    #endif
        }

    }
    //  G4cerr <<"Finished capture loop "<<G4endl;

    //G4cerr <<"pre- DoTimeStep 4"<<G4endl;
    DoTimeStep(DBL_MAX);
    //G4cerr <<"post- DoTimeStep 4"<<G4endl;


}

//----------------------------------------------------------------------------
G4double G4BinaryCascade::CorrectShortlivedPrimaryForFermi(
        G4KineticTrack* primary,G4KineticTrackVector target_collection)
//----------------------------------------------------------------------------
{
    G4double Efermi(0);
    if (primary->GetState() == G4KineticTrack::inside ) {
        G4int PDGcode=primary->GetDefinition()->GetPDGEncoding();
        Efermi=((G4RKPropagation *)thePropagator)->GetField(PDGcode,primary->GetPosition());

        if ( std::abs(PDGcode) > 1000 && PDGcode != 2112 && PDGcode != 2212 )
        {
            Efermi = ((G4RKPropagation *)thePropagator)->GetField(G4Neutron::Neutron()->GetPDGEncoding(),primary->GetPosition());
            G4LorentzVector mom4Primary=primary->Get4Momentum();
            primary->Update4Momentum(mom4Primary.e() - Efermi);
        }

        std::vector<G4KineticTrack *>::iterator titer;
        for ( titer=target_collection.begin() ; titer!=target_collection.end(); ++titer)
        {
            const G4ParticleDefinition * aDef=(*titer)->GetDefinition();
            G4int aCode=aDef->GetPDGEncoding();
            G4ThreeVector aPos=(*titer)->GetPosition();
            Efermi+= ((G4RKPropagation *)thePropagator)->GetField(aCode, aPos);
        }
    }
    return Efermi;
}

//----------------------------------------------------------------------------
G4bool G4BinaryCascade::CorrectShortlivedFinalsForFermi(G4KineticTrackVector * products,
        G4double initial_Efermi)
//----------------------------------------------------------------------------
{
    G4double final_Efermi(0);
    G4KineticTrackVector resonances;
    for ( std::vector<G4KineticTrack *>::iterator i =products->begin(); i != products->end(); i++)
    {
        G4int PDGcode=(*i)->GetDefinition()->GetPDGEncoding();
        //       G4cout << " PDGcode, state " << PDGcode << " " << (*i)->GetState()<<G4endl;
        final_Efermi+=((G4RKPropagation *)thePropagator)->GetField(PDGcode,(*i)->GetPosition());
        if ( std::abs(PDGcode) > 1000 && PDGcode != 2112 && PDGcode != 2212 )
        {
            resonances.push_back(*i);
        }
    }
    if ( resonances.size() > 0 )
    {
        G4double delta_Fermi= (initial_Efermi-final_Efermi)/resonances.size();
        for (std::vector<G4KineticTrack *>::iterator res=resonances.begin(); res != resonances.end(); res++)
        {
            G4LorentzVector mom=(*res)->Get4Momentum();
            G4double mass2=mom.mag2();
            G4double newEnergy=mom.e() + delta_Fermi;
            G4double newEnergy2= newEnergy*newEnergy;
            //G4cout << "mom = " << mom <<" newE " << newEnergy<< G4endl;
            if ( newEnergy2 < mass2 )
            {
                return false;
            }
            G4ThreeVector mom3=std::sqrt(newEnergy2 - mass2) * mom.vect().unit();
            (*res)->Set4Momentum(G4LorentzVector(mom3,newEnergy));
            	  //G4cout << " correct resonance from /to " << mom.e() << " / " << newEnergy<<
            	  //		    " 3mom from/to " << mom.vect() << " / " << mom3 << G4endl;
        }
    }
    return true;
}

//----------------------------------------------------------------------------
void G4BinaryCascade::CorrectFinalPandE()
//----------------------------------------------------------------------------
//
//  Modify momenta of outgoing particles.
//   Assume two body decay, nucleus(@nominal mass) + sum of final state particles(SFSP).
//   momentum of SFSP shall be less than momentum for two body decay.
//
{
#ifdef debug_BIC_CorrectFinalPandE
    G4cerr << "BIC: -CorrectFinalPandE called" << G4endl;
#endif

    if ( theFinalState.size() == 0 ) return;

    G4KineticTrackVector::iterator i;
    G4LorentzVector pNucleus=GetFinal4Momentum();
    if ( pNucleus.e() == 0 ) return;    // check against explicit 0 from GetNucleus4Momentum()
#ifdef debug_BIC_CorrectFinalPandE
    G4cerr << " -CorrectFinalPandE 3" << G4endl;
#endif
    G4LorentzVector pFinals(0);
    G4int nFinals(0);
    for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
    {
        pFinals += (*i)->Get4Momentum();
        ++nFinals;
#ifdef debug_BIC_CorrectFinalPandE
        G4cout <<"CorrectFinalPandE a final " << (*i)->GetDefinition()->GetParticleName()
        		           << " 4mom " << (*i)->Get4Momentum()<< G4endl;
#endif
    }
#ifdef debug_BIC_CorrectFinalPandE
    G4cout << "CorrectFinalPandE pN pF: " <<pNucleus << " " <<pFinals << G4endl;
#endif
    G4LorentzVector pCM=pNucleus + pFinals;

    G4LorentzRotation toCMS(-pCM.boostVector());
    pFinals *=toCMS;
#ifdef debug_BIC_CorrectFinalPandE
    G4cout << "CorrectFinalPandE pCM, CMS pCM " << pCM << " " <<toCMS*pCM<< G4endl;
    G4cout << "CorrectFinal CMS pN pF " <<toCMS*pNucleus << " "
            <<pFinals << G4endl
            << " nucleus initial mass : " <<GetFinal4Momentum().mag()
            <<" massInNucleus m(nucleus) m(finals) std::sqrt(s): " << massInNucleus << " " <<pNucleus.mag()<< " "
            << pFinals.mag() << " " << pCM.mag() << G4endl;
#endif

    G4LorentzRotation toLab = toCMS.inverse();

    G4double s0 = pCM.mag2();
    G4double m10 = GetIonMass(currentZ,currentA);
    G4double m20 = pFinals.mag();
    if( s0-(m10+m20)*(m10+m20) < 0 )
    {
#ifdef debug_BIC_CorrectFinalPandE
        G4cout << "G4BinaryCascade::CorrectFinalPandE() : error! " << G4endl;

        G4cout << "not enough mass to correct: mass^2, A,Z, mass(nucl), mass(finals) "
                << (s0-(m10+m20)*(m10+m20)) << " "
                << currentA << " " << currentZ << " "
                << m10 << " " << m20
                << G4endl;
        G4cerr << " -CorrectFinalPandE 4" << G4endl;

        PrintKTVector(&theFinalState," mass problem");
#endif
        return;
    }

    // Three momentum in cm system
    G4double pInCM = std::sqrt((s0-(m10+m20)*(m10+m20))*(s0-(m10-m20)*(m10-m20))/(4.*s0));
#ifdef debug_BIC_CorrectFinalPandE
    G4cout <<" CorrectFinalPandE pInCM  new, CURRENT, ratio : " << pInCM
            << " " << (pFinals).vect().mag()<< " " <<  pInCM/(pFinals).vect().mag() << G4endl;
#endif
    if ( pFinals.vect().mag() > pInCM )
    {
        G4ThreeVector p3finals=pInCM*pFinals.vect().unit();

        //    G4ThreeVector deltap=(p3finals - pFinals.vect() ) / nFinals;
        G4double factor=std::max(0.98,pInCM/pFinals.vect().mag());   // small correction
        G4LorentzVector qFinals(0);
        for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
        {
            //      G4ThreeVector p3((toCMS*(*i)->Get4Momentum()).vect() + deltap);
            G4ThreeVector p3(factor*(toCMS*(*i)->Get4Momentum()).vect());
            G4LorentzVector p(p3,std::sqrt((*i)->Get4Momentum().mag2() + p3.mag2()));
            qFinals += p;
            p *= toLab;
#ifdef debug_BIC_CorrectFinalPandE
            G4cout << " final p corrected: " << p << G4endl;
#endif
            (*i)->Set4Momentum(p);
        }
#ifdef debug_BIC_CorrectFinalPandE
        G4cout << "CorrectFinalPandE nucleus corrected mass : " << GetFinal4Momentum() << " "
                <<GetFinal4Momentum().mag() << G4endl
                << " CMS pFinals , mag, 3.mag : " << qFinals << " " << qFinals.mag() << " " << qFinals.vect().mag()<< G4endl;
        G4cerr << " -CorrectFinalPandE 5 " << factor <<  G4endl;
#endif
    }
#ifdef debug_BIC_CorrectFinalPandE
    else { G4cerr << " -CorrectFinalPandE 6 - no correction done" << G4endl; }
#endif

}

//----------------------------------------------------------------------------
void G4BinaryCascade::UpdateTracksAndCollisions(
        //----------------------------------------------------------------------------
        G4KineticTrackVector * oldSecondaries,
        G4KineticTrackVector * oldTarget,
        G4KineticTrackVector * newSecondaries)
{
    std::vector<G4KineticTrack *>::iterator iter1, iter2;

    // remove old secondaries from the secondary list
    if(oldSecondaries)
    {
        if(!oldSecondaries->empty())
        {
            for(iter1 = oldSecondaries->begin(); iter1 != oldSecondaries->end();
                    ++iter1)
            {
                iter2 = std::find(theSecondaryList.begin(), theSecondaryList.end(),
                        *iter1);
                if ( iter2 != theSecondaryList.end() ) theSecondaryList.erase(iter2);
            }
            theCollisionMgr->RemoveTracksCollisions(oldSecondaries);
        }
    }

    // remove old target from the target list
    if(oldTarget)
    {
        // G4cout << "################## Debugging 0 "<<G4endl;
        if(oldTarget->size()!=0)
        {

            // G4cout << "################## Debugging 1 "<<oldTarget->size()<<G4endl;
            for(iter1 = oldTarget->begin(); iter1 != oldTarget->end(); ++iter1)
            {
                iter2 = std::find(theTargetList.begin(), theTargetList.end(),
                        *iter1);
                theTargetList.erase(iter2);
            }
            theCollisionMgr->RemoveTracksCollisions(oldTarget);
        }
    }

    if(newSecondaries)
    {
        if(!newSecondaries->empty())
        {
            // insert new secondaries in the secondary list
            for(iter1 = newSecondaries->begin(); iter1 != newSecondaries->end();
                    ++iter1)
            {
                theSecondaryList.push_back(*iter1);
                if ((*iter1)->GetState() == G4KineticTrack::undefined)
                {
                   PrintKTVector(*iter1, "undefined in FindCollisions");
                }


            }
            // look for collisions of new secondaries
            FindCollisions(newSecondaries);
        }
    }
    // G4cout << "Exiting ... "<<oldTarget<<G4endl;
}


class SelectFromKTV
{
private:
    G4KineticTrackVector * ktv;
    G4KineticTrack::CascadeState wanted_state;
public:
    SelectFromKTV(G4KineticTrackVector * out, G4KineticTrack::CascadeState astate)
    :
        ktv(out), wanted_state(astate)
    {};
    void operator() (G4KineticTrack *& kt) const
    {
        if ( (kt)->GetState() == wanted_state ) ktv->push_back(kt);
    };
};



//----------------------------------------------------------------------------
G4bool G4BinaryCascade::DoTimeStep(G4double theTimeStep)
//----------------------------------------------------------------------------
{

#ifdef debug_BIC_DoTimeStep
    G4ping debug("debug_G4BinaryCascade");
    debug.push_back("======> DoTimeStep 1"); debug.dump();
    G4cerr <<"G4BinaryCascade::DoTimeStep: enter step="<< theTimeStep
            << " , time="<<theCurrentTime << G4endl;
    PrintKTVector(&theSecondaryList, std::string("DoTimeStep - theSecondaryList"));
    //PrintKTVector(&theTargetList, std::string("DoTimeStep - theTargetList"));
#endif

    G4bool success=true;
    std::vector<G4KineticTrack *>::iterator iter;

    G4KineticTrackVector * kt_outside = new G4KineticTrackVector;
    std::for_each( theSecondaryList.begin(),theSecondaryList.end(),
            SelectFromKTV(kt_outside,G4KineticTrack::outside));
    //PrintKTVector(kt_outside, std::string("DoTimeStep - found outside"));

    G4KineticTrackVector * kt_inside = new G4KineticTrackVector;
    std::for_each( theSecondaryList.begin(),theSecondaryList.end(),
            SelectFromKTV(kt_inside, G4KineticTrack::inside));
    //  PrintKTVector(kt_inside, std::string("DoTimeStep - found inside"));
    //-----
    G4KineticTrackVector dummy;   // needed for re-usability
#ifdef debug_BIC_DoTimeStep
    G4cout << "NOW WE ARE ENTERING THE TRANSPORT"<<G4endl;
#endif

    // =================== Here we move the particles  ===================

    thePropagator->Transport(theSecondaryList, dummy, theTimeStep);

    // =================== Here we move the particles  ===================

    //------

    theMomentumTransfer += thePropagator->GetMomentumTransfer();
#ifdef debug_BIC_DoTimeStep
    G4cout << "DoTimeStep : theMomentumTransfer = " << theMomentumTransfer << G4endl;
    PrintKTVector(&theSecondaryList, std::string("DoTimeStep - secondaries aft trsprt"));
#endif

    //_DebugEpConservation(" after stepping");

    // Partclies which went INTO nucleus

    G4KineticTrackVector * kt_gone_in = new G4KineticTrackVector;
    std::for_each( kt_outside->begin(),kt_outside->end(),
            SelectFromKTV(kt_gone_in,G4KineticTrack::inside));
    //  PrintKTVector(kt_gone_in, std::string("DoTimeStep - gone in"));


    // Partclies which  went OUT OF nucleus
    G4KineticTrackVector * kt_gone_out = new G4KineticTrackVector;
    std::for_each( kt_inside->begin(),kt_inside->end(),
            SelectFromKTV(kt_gone_out, G4KineticTrack::gone_out));

    //  PrintKTVector(kt_gone_out, std::string("DoTimeStep - gone out"));

    G4KineticTrackVector *fail=CorrectBarionsOnBoundary(kt_gone_in,kt_gone_out);

    if ( fail )
    {
        // some particle(s) supposed to enter/leave were miss_nucleus/captured by the correction
        kt_gone_in->clear();
        std::for_each( kt_outside->begin(),kt_outside->end(),
                SelectFromKTV(kt_gone_in,G4KineticTrack::inside));

        kt_gone_out->clear();
        std::for_each( kt_inside->begin(),kt_inside->end(),
                SelectFromKTV(kt_gone_out, G4KineticTrack::gone_out));

#ifdef debug_BIC_DoTimeStep
        PrintKTVector(fail,std::string(" Failed to go in/out -> miss_nucleus/captured"));
        PrintKTVector(kt_gone_in, std::string("recreated kt_gone_in"));
        PrintKTVector(kt_gone_out, std::string("recreated kt_gone_out"));
#endif
        delete fail;
    }

    // Add tracks missing nucleus and tracks going straight though  to addFinals
    std::for_each( kt_outside->begin(),kt_outside->end(),
            SelectFromKTV(kt_gone_out,G4KineticTrack::miss_nucleus));
    //PrintKTVector(kt_gone_out, std::string("miss to append to final state.."));
    std::for_each( kt_outside->begin(),kt_outside->end(),
            SelectFromKTV(kt_gone_out,G4KineticTrack::gone_out));

#ifdef debug_BIC_DoTimeStep
    PrintKTVector(kt_gone_out, std::string("append gone_outs to final state.. theFinalState"));
#endif

    theFinalState.insert(theFinalState.end(),
            kt_gone_out->begin(),kt_gone_out->end());

    // Partclies which could not leave nucleus,  captured...
    G4KineticTrackVector * kt_captured = new G4KineticTrackVector;
    std::for_each( theSecondaryList.begin(),theSecondaryList.end(),
            SelectFromKTV(kt_captured, G4KineticTrack::captured));

    // Check no track is part in next collision, ie.
    //  this step was to far, and collisions should not occur any more

    if ( theCollisionMgr->Entries()> 0 )
    {
        if (kt_gone_out->size() )
        {
            G4KineticTrack * nextPrimary = theCollisionMgr->GetNextCollision()->GetPrimary();
            iter = std::find(kt_gone_out->begin(),kt_gone_out->end(),nextPrimary);
            if ( iter !=  kt_gone_out->end() )
            {
                success=false;
#ifdef debug_BIC_DoTimeStep
                G4cout << " DoTimeStep - WARNING: deleting current collision!" << G4endl;
#endif
            }
        }
        if ( kt_captured->size() )
        {
            G4KineticTrack * nextPrimary = theCollisionMgr->GetNextCollision()->GetPrimary();
            iter = std::find(kt_captured->begin(),kt_captured->end(),nextPrimary);
            if ( iter !=  kt_captured->end() )
            {
                success=false;
#ifdef debug_BIC_DoTimeStep
                G4cout << " DoTimeStep - WARNING: deleting current collision!" << G4endl;
#endif
            }
        }

    }
    // PrintKTVector(kt_gone_out," kt_gone_out be4 updatetrack...");
    UpdateTracksAndCollisions(kt_gone_out,0 ,0);


    if ( kt_captured->size() )
    {
        theCapturedList.insert(theCapturedList.end(),
                kt_captured->begin(),kt_captured->end());
        //should be      std::for_each(kt_captured->begin(),kt_captured->end(),
        //              std::mem_fun(&G4KineticTrack::Hit));
        // but VC 6 requires:
        std::vector<G4KineticTrack *>::iterator i_captured;
        for(i_captured=kt_captured->begin();i_captured!=kt_captured->end();i_captured++)
        {
            (*i_captured)->Hit();
        }
        //     PrintKTVector(kt_captured," kt_captured be4 updatetrack...");
        UpdateTracksAndCollisions(kt_captured, NULL, NULL);
    }

#ifdef debug_G4BinaryCascade
    delete kt_inside;
    kt_inside = new G4KineticTrackVector;
    std::for_each( theSecondaryList.begin(),theSecondaryList.end(),
            SelectFromKTV(kt_inside, G4KineticTrack::inside));
    if ( currentZ != (GetTotalCharge(theTargetList)
            + GetTotalCharge(theCapturedList)
            + GetTotalCharge(*kt_inside)) )
    {
        G4cout << " error-DoTimeStep aft, A, Z: " << currentA << " " << currentZ
                << " sum(tgt,capt,active) "
                << GetTotalCharge(theTargetList) + GetTotalCharge(theCapturedList) + GetTotalCharge(*kt_inside)
                << " targets: "  << GetTotalCharge(theTargetList)
                << " captured: " << GetTotalCharge(theCapturedList)
                << " active: "   << GetTotalCharge(*kt_inside)
                << G4endl;
    }
#endif

    delete kt_inside;
    delete kt_outside;
    delete kt_captured;
    delete kt_gone_in;
    delete kt_gone_out;

    //  G4cerr <<"G4BinaryCascade::DoTimeStep: exit "<<G4endl;
    theCurrentTime += theTimeStep;

    //debug.push_back("======> DoTimeStep 2"); debug.dump();
    return success;

}

//----------------------------------------------------------------------------
G4KineticTrackVector* G4BinaryCascade::CorrectBarionsOnBoundary(
        G4KineticTrackVector *in,
        G4KineticTrackVector *out)
//----------------------------------------------------------------------------
{
    G4KineticTrackVector * kt_fail(0);
    std::vector<G4KineticTrack *>::iterator iter;
    //  G4cout << "CorrectBarionsOnBoundary,currentZ,currentA,"
    //         << currentZ << " "<< currentA << G4endl;
    if (in->size())
    {
        G4int secondaries_in(0);
        G4int secondaryBarions_in(0);
        G4int secondaryCharge_in(0);
        G4double secondaryMass_in(0);

        for ( iter =in->begin(); iter != in->end(); ++iter)
        {
            ++secondaries_in;
            secondaryCharge_in += G4lrint((*iter)->GetDefinition()->GetPDGCharge()/eplus);
            if ((*iter)->GetDefinition()->GetBaryonNumber()!=0 )
            {
                secondaryBarions_in += (*iter)->GetDefinition()->GetBaryonNumber();
                if((*iter)->GetDefinition() == G4Neutron::Neutron() ||
                        (*iter)->GetDefinition() == G4Proton::Proton()  )
                {
                    secondaryMass_in += (*iter)->GetDefinition()->GetPDGMass();
                } else 	  {
                    secondaryMass_in += G4Proton::Proton()->GetPDGMass();
                }
            }
        }
        G4double mass_initial= GetIonMass(currentZ,currentA);

        currentZ += secondaryCharge_in;
        currentA += secondaryBarions_in;

        //  G4cout << "CorrectBarionsOnBoundary,secondaryCharge_in, secondaryBarions_in "
        //         <<    secondaryCharge_in << " "<<  secondaryBarions_in << G4endl;

        G4double mass_final= GetIonMass(currentZ,currentA);

        G4double correction= secondaryMass_in + mass_initial - mass_final;
        if (secondaries_in>1)
        {correction /= secondaries_in;}

#ifdef debug_BIC_CorrectBarionsOnBoundary
        G4cout << "CorrectBarionsOnBoundary,currentZ,currentA,"
                << "secondaryCharge_in,secondaryBarions_in,"
                << "energy correction,m_secondry,m_nucl_init,m_nucl_final "
                << currentZ << " "<< currentA <<" "
                << secondaryCharge_in<<" "<<secondaryBarions_in<<" "
                << correction << " "
                << secondaryMass_in << " "
                << mass_initial << " "
                << mass_final << " "
                << G4endl;
        PrintKTVector(in,std::string("in be4 correction"));
#endif
        for ( iter = in->begin(); iter != in->end(); ++iter)
        {
            if ((*iter)->GetTrackingMomentum().e()+correction > (*iter)->GetActualMass())
            {
                (*iter)->UpdateTrackingMomentum((*iter)->GetTrackingMomentum().e() + correction);
            } else {
                //particle cannot go in, put to miss_nucleus
                G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;
                (*iter)->SetState(G4KineticTrack::miss_nucleus);
                // Undo correction for Colomb Barrier
                G4double barrier=RKprop->GetBarrier((*iter)->GetDefinition()->GetPDGEncoding());
                (*iter)->UpdateTrackingMomentum((*iter)->GetTrackingMomentum().e() + barrier);
                if ( ! kt_fail ) kt_fail=new G4KineticTrackVector;
                kt_fail->push_back(*iter);
                currentZ -= G4lrint((*iter)->GetDefinition()->GetPDGCharge()/eplus);
                currentA -= (*iter)->GetDefinition()->GetBaryonNumber();

            }

        }
#ifdef debug_BIC_CorrectBarionsOnBoundary
        G4cout << " CorrectBarionsOnBoundary, aft, Z, A, sec-Z,A,m,m_in_nucleus "
                << currentZ << " " << currentA << " "
                << secondaryCharge_in << " " << secondaryBarions_in << " "
                << secondaryMass_in  << " "
                << G4endl;
        PrintKTVector(in,std::string("in AFT correction"));
#endif

    }
    //----------------------------------------------
    if (out->size())
    {
        G4int secondaries_out(0);
        G4int secondaryBarions_out(0);
        G4int secondaryCharge_out(0);
        G4double secondaryMass_out(0);

        for ( iter =out->begin(); iter != out->end(); ++iter)
        {
            ++secondaries_out;
            secondaryCharge_out += G4lrint((*iter)->GetDefinition()->GetPDGCharge()/eplus);
            if ((*iter)->GetDefinition()->GetBaryonNumber() !=0 )
            {
                secondaryBarions_out += (*iter)->GetDefinition()->GetBaryonNumber();
                if((*iter)->GetDefinition() == G4Neutron::Neutron() ||
                        (*iter)->GetDefinition() == G4Proton::Proton()  )
                {
                    secondaryMass_out += (*iter)->GetDefinition()->GetPDGMass();
                } else {
                    secondaryMass_out += G4Neutron::Neutron()->GetPDGMass();
                }
            }
        }

        G4double mass_initial=  GetIonMass(currentZ,currentA);
        currentA -=secondaryBarions_out;
        currentZ -=secondaryCharge_out;

        //  G4cout << "CorrectBarionsOnBoundary,secondaryCharge_out, secondaryBarions_out "
        //         <<    secondaryCharge_out << " "<<  secondaryBarions_out << G4endl;

        //                        a delta minus will do currentZ < 0 in light nuclei
        //     if (currentA < 0 || currentZ < 0 )
        if (currentA < 0 )
        {
            G4cerr << "G4BinaryCascade - secondaryBarions_out,secondaryCharge_out " <<
                    secondaryBarions_out << " " << secondaryCharge_out << G4endl;
            PrintKTVector(&theTargetList,"CorrectBarionsOnBoundary Target");
            PrintKTVector(&theCapturedList,"CorrectBarionsOnBoundary Captured");
            PrintKTVector(&theSecondaryList,"CorrectBarionsOnBoundary Secondaries");
            G4cerr << "G4BinaryCascade - currentA, currentZ " << currentA << " " << currentZ << G4endl;
            throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::CorrectBarionsOnBoundary() - fatal error");
        }
        G4double mass_final=GetIonMass(currentZ,currentA);
        G4double correction= mass_initial - mass_final - secondaryMass_out;
        // G4cout << "G4BinaryCascade::CorrectBarionsOnBoundary() total out correction: " << correction << G4endl;

        if (secondaries_out>1) correction /= secondaries_out;
#ifdef debug_BIC_CorrectBarionsOnBoundary
        G4cout << "DoTimeStep,(current Z,A),"
                << "(secondaries out,Charge,Barions),"
                <<"* energy correction,(m_secondry,m_nucl_init,m_nucl_final) "
                << "("<< currentZ << ","<< currentA <<") ("
                << secondaries_out << ","
                << secondaryCharge_out<<","<<secondaryBarions_out<<") * "
                << correction << " ("
                << secondaryMass_out << ", "
                << mass_initial << ", "
                << mass_final << ")"
                << G4endl;
        PrintKTVector(out,std::string("out be4 correction"));
#endif

        for ( iter = out->begin(); iter != out->end(); ++iter)
        {
            if ((*iter)->GetTrackingMomentum().e()+correction > (*iter)->GetActualMass())
            {
                (*iter)->UpdateTrackingMomentum((*iter)->GetTrackingMomentum().e() + correction);
            } else
            {
                // particle cannot go out due to change of nuclear potential!
                //  capture protons and neutrons;
                if(((*iter)->GetDefinition() == G4Proton::Proton()) ||
                        ((*iter)->GetDefinition() == G4Neutron::Neutron()))
                {
                    (*iter)->SetState(G4KineticTrack::captured);
                    // Undo correction for Colomb Barrier
                    G4double barrier=((G4RKPropagation *)thePropagator)->GetBarrier((*iter)->GetDefinition()->GetPDGEncoding());
                    (*iter)->UpdateTrackingMomentum((*iter)->GetTrackingMomentum().e() - barrier);
                    if ( kt_fail == 0 ) kt_fail=new G4KineticTrackVector;
                    kt_fail->push_back(*iter);
                    currentZ += G4lrint((*iter)->GetDefinition()->GetPDGCharge()/eplus);
                    currentA += (*iter)->GetDefinition()->GetBaryonNumber();
                }
#ifdef debug_BIC_CorrectBarionsOnBoundary
                else
                {
                    G4cout << "Not correcting outgoing " << *iter << " "
                            << (*iter)->GetDefinition()->GetPDGEncoding() << " "
                            << (*iter)->GetDefinition()->GetParticleName() << G4endl;
                    PrintKTVector(out,std::string("outgoing, one not corrected"));
                }
#endif
            }
        }

#ifdef debug_BIC_CorrectBarionsOnBoundary
        PrintKTVector(out,std::string("out AFTER correction"));
        G4cout << " DoTimeStep, nucl-update, A, Z, sec-Z,A,m,m_in_nucleus, table-mass, delta "
                << currentA << " "<< currentZ << " "
                << secondaryCharge_out << " "<< secondaryBarions_out << " "<<
                secondaryMass_out << " "
                << massInNucleus << " "
                << GetIonMass(currentZ,currentA)
                << " " << massInNucleus - GetIonMass(currentZ,currentA)
                << G4endl;
#endif
    }

    return kt_fail;
}


//----------------------------------------------------------------------------

G4Fragment * G4BinaryCascade::FindFragments()
//----------------------------------------------------------------------------
{

#ifdef debug_BIC_FindFragments
    G4cout << "target, captured, secondary: "
            << theTargetList.size() << " "
            << theCapturedList.size()<< " "
            << theSecondaryList.size()
            << G4endl;
#endif

    G4int a = theTargetList.size()+theCapturedList.size();
    G4int zTarget = 0;
    G4KineticTrackVector::iterator i;
    for(i = theTargetList.begin(); i != theTargetList.end(); ++i)
    {
        if(G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus) == 1 )
        {
            zTarget++;
        }
    }

    G4int zCaptured = 0;
    G4LorentzVector CapturedMomentum(0.,0.,0.,0.);
    for(i = theCapturedList.begin(); i != theCapturedList.end(); ++i)
    {
        CapturedMomentum += (*i)->Get4Momentum();
        if(G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus) == 1 )
        {
            zCaptured++;
        }
    }

    G4int z = zTarget+zCaptured;

#ifdef debug_G4BinaryCascade
    if ( z != (GetTotalCharge(theTargetList) + GetTotalCharge(theCapturedList)) )
    {
        G4cout << " FindFragment Counting error z a " << z << " " <<a << " "
                << GetTotalCharge(theTargetList) << " " <<  GetTotalCharge(theCapturedList)<<
                G4endl;
        PrintKTVector(&theTargetList, std::string("theTargetList"));
        PrintKTVector(&theCapturedList, std::string("theCapturedList"));
    }
#endif
    //debug
    /*
     *   G4cout << " Fragment mass table / real "
     *          << GetIonMass(z, a)
     * 	 << " / " << GetFinal4Momentum().mag()
     * 	 << " difference "
     * 	 <<  GetFinal4Momentum().mag() - GetIonMass(z, a)
     * 	 << G4endl;
     */
    //
    //  if(getenv("BCDEBUG") ) G4cerr << "Fragment A, Z "<< a <<" "<< z<<G4endl;
    if ( z < 1 ) return 0;

    G4int holes = the3DNucleus->GetMassNumber() - theTargetList.size();
    G4int excitons = theCapturedList.size();
#ifdef debug_BIC_FindFragments
    G4cout << "Fragment: a= " << a << " z= " << z << " particles= " <<  excitons
            << " Charged= " << zCaptured << " holes= " << holes
            << " excitE= " <<GetExcitationEnergy()
            << " Final4Momentum= " << GetFinalNucleusMomentum() << " capturMomentum= " << CapturedMomentum
            << G4endl;
#endif

    G4Fragment * fragment = new G4Fragment(a,z,GetFinalNucleusMomentum());
    fragment->SetNumberOfHoles(holes);

    //GF  fragment->SetNumberOfParticles(excitons-holes);
    fragment->SetNumberOfParticles(excitons);
    fragment->SetNumberOfCharged(zCaptured);

    return fragment;
}

//----------------------------------------------------------------------------

G4LorentzVector G4BinaryCascade::GetFinal4Momentum()
//----------------------------------------------------------------------------
// Return momentum of reminder nulceus;
//  ie. difference of (initial state(primaries+nucleus) - final state) particles, ignoring remnant nucleus
{
    G4LorentzVector final4Momentum = theInitial4Mom + theProjectile4Momentum;
    G4LorentzVector finals(0,0,0,0);
    for(G4KineticTrackVector::iterator i = theFinalState.begin(); i != theFinalState.end(); ++i)
    {
        final4Momentum -= (*i)->Get4Momentum();
        finals		   += (*i)->Get4Momentum();
    }

    if(final4Momentum.e()> 0 && (final4Momentum.vect()/final4Momentum.e()).mag()>1.0 && currentA > 0)
    {
#ifdef debug_BIC_Final4Momentum
        G4cerr << G4endl;
        G4cerr << "G4BinaryCascade::GetFinal4Momentum - Fatal"<<G4endl;
        G4KineticTrackVector::iterator i;
        G4cerr <<"Total initial 4-momentum " << theProjectile4Momentum << G4endl;
        G4cerr <<" GetFinal4Momentum: Initial nucleus "<<theInitial4Mom<<G4endl;
        for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
        {
            G4cerr <<" Final state: "<<(*i)->Get4Momentum()<<(*i)->GetDefinition()->GetParticleName()<<G4endl;
        }
        G4cerr << "Sum( 4-mom ) finals " << finals << G4endl;
        G4cerr<< " Final4Momentum = "<<final4Momentum <<" "<<final4Momentum.m()<<G4endl;
        G4cerr <<" current A, Z = "<< currentA<<", "<<currentZ<<G4endl;
        G4cerr << G4endl;
#endif

        final4Momentum=G4LorentzVector(0,0,0,0);
    }
    return final4Momentum;
}

//----------------------------------------------------------------------------
G4LorentzVector G4BinaryCascade::GetFinalNucleusMomentum()
//----------------------------------------------------------------------------
{
    // return momentum of nucleus for use with precompound model; also keep transformation to
    //   apply to precompoud products.

    G4LorentzVector CapturedMomentum(0,0,0,0);
    G4KineticTrackVector::iterator i;
    //  G4cout << "GetFinalNucleusMomentum Captured size: " <<theCapturedList.size() << G4endl;
    for(i = theCapturedList.begin(); i != theCapturedList.end(); ++i)
    {
        CapturedMomentum += (*i)->Get4Momentum();
    }
    //G4cout << "GetFinalNucleusMomentum CapturedMomentum= " <<CapturedMomentum << G4endl;
    //  G4cerr << "it 9"<<G4endl;

    G4LorentzVector NucleusMomentum = GetFinal4Momentum();
    if ( NucleusMomentum.e() > 0 )
    {
        // G4cout << "GetFinalNucleusMomentum GetFinal4Momentum= " <<NucleusMomentum <<" "<<NucleusMomentum.mag()<<G4endl;
        // boost nucleus to a frame such that the momentum of nucleus == momentum of Captured
        G4ThreeVector boost= (NucleusMomentum.vect() -CapturedMomentum.vect())/NucleusMomentum.e();
        if(boost.mag2()>1.0)
        {
#     ifdef debug_BIC_FinalNucleusMomentum
            G4cerr << "G4BinaryCascade::GetFinalNucleusMomentum - Fatal"<<G4endl;
            G4cerr << "it 0"<<boost <<G4endl;
            G4cerr << "it 01"<<NucleusMomentum<<" "<<CapturedMomentum<<" "<<G4endl;
            G4cout <<" testing boost "<<boost<<" "<<boost.mag()<<G4endl;
#      endif
            boost=G4ThreeVector(0,0,0);
            NucleusMomentum=G4LorentzVector(0,0,0,0);
        }
        G4LorentzRotation  nucleusBoost( -boost );
        precompoundLorentzboost.set( boost );
#ifdef debug_debug_BIC_FinalNucleusMomentum
        G4cout << "GetFinalNucleusMomentum be4 boostNucleusMomentum, CapturedMomentum"<<NucleusMomentum<<" "<<CapturedMomentum<<" "<<G4endl;
#endif
        NucleusMomentum *= nucleusBoost;
#ifdef debug_BIC_FinalNucleusMomentum
        G4cout << "GetFinalNucleusMomentum aft boost GetFinal4Momentum= " <<NucleusMomentum <<G4endl;
#endif
    }
    return NucleusMomentum;
}

//----------------------------------------------------------------------------
G4ReactionProductVector * G4BinaryCascade::Propagate1H1(
        //----------------------------------------------------------------------------
        G4KineticTrackVector * secondaries, G4V3DNucleus * nucleus)
{
    G4ReactionProductVector * products = new G4ReactionProductVector;
    const G4ParticleDefinition * aHTarg = G4Proton::ProtonDefinition();
    if (nucleus->GetCharge() == 0) aHTarg = G4Neutron::NeutronDefinition();
    G4double mass = aHTarg->GetPDGMass();
    G4KineticTrackVector * secs = 0;
    G4ThreeVector pos(0,0,0);
    G4LorentzVector mom(mass);
    G4KineticTrack aTarget(aHTarg, 0., pos, mom);
    G4bool done(false);
    std::vector<G4KineticTrack *>::iterator iter, jter;
    // data member    static G4Scatterer theH1Scatterer;
    //G4cout << " start 1H1 for " << (*secondaries).front()->GetDefinition()->GetParticleName()
    //       << " on " << aHTarg->GetParticleName() << G4endl;
    G4int tryCount(0);
    while(!done && tryCount++ <200)                                 /* Loop checking, 31.08.2015, G.Folger */
    {
        if(secs)
        {
            std::for_each(secs->begin(), secs->end(), DeleteKineticTrack());
            delete secs;
        }
        secs = theH1Scatterer->Scatter(*(*secondaries).front(), aTarget);
#ifdef debug_H1_BinaryCascade
        PrintKTVector(secs," From Scatter");
#endif
        for(size_t ss=0; secs && ss<secs->size(); ss++)
        {
            // must have one resonance in final state, or it was elastic, not allowed here.
            if((*secs)[ss]->GetDefinition()->IsShortLived()) done = true;
        }
    }

    ClearAndDestroy(&theFinalState);
    ClearAndDestroy(secondaries);
    delete secondaries;

    for(size_t current=0; secs && current<secs->size(); current++)
    {
        if((*secs)[current]->GetDefinition()->IsShortLived())
        {
            done = true; 	// must have one resonance in final state, elastic not allowed here!
            G4KineticTrackVector * dec = (*secs)[current]->Decay();
            for(jter=dec->begin(); jter != dec->end(); jter++)
            {
                //G4cout << "Decay"<<G4endl;
                secs->push_back(*jter);
                //G4cout << "decay "<<(*jter)->GetDefinition()->GetParticleName()<<G4endl;
            }
            delete (*secs)[current];
            delete dec;
        }
        else
        {
            //G4cout << "FS "<<G4endl;
            //G4cout << "FS "<<(*secs)[current]->GetDefinition()->GetParticleName()<<G4endl;
            theFinalState.push_back((*secs)[current]);
        }
    }

    delete secs;
#ifdef debug_H1_BinaryCascade
    PrintKTVector(&theFinalState," FinalState");
#endif
    for(iter = theFinalState.begin(); iter != theFinalState.end(); ++iter)
    {
        G4KineticTrack * kt = *iter;
        G4ReactionProduct * aNew = new G4ReactionProduct(kt->GetDefinition());
        aNew->SetMomentum(kt->Get4Momentum().vect());
        aNew->SetTotalEnergy(kt->Get4Momentum().e());
        aNew->SetCreatorModel(theBIC_ID);
        products->push_back(aNew);
#ifdef debug_H1_BinaryCascade
        if (! kt->GetDefinition()->GetPDGStable() )
        {
            if (kt->GetDefinition()->IsShortLived())
            {
                G4cout << "final shortlived : ";
            } else
            {
                G4cout << "final un stable : ";
            }
            G4cout <<kt->GetDefinition()->GetParticleName()<< G4endl;
        }
#endif
        delete kt;
    }
    theFinalState.clear();
    return products;

}

//----------------------------------------------------------------------------
G4ThreeVector G4BinaryCascade::GetSpherePoint(
        G4double r, const G4LorentzVector & mom4)
//----------------------------------------------------------------------------
{
    //  Get a point outside radius.
    //     point is random in plane (circle of radius r) orthogonal to mom,
    //      plus -1*r*mom->vect()->unit();
    G4ThreeVector o1, o2;
    G4ThreeVector mom = mom4.vect();

    o1= mom.orthogonal();  // we simply need any vector non parallel
    o2= mom.cross(o1);     //  o2 is now orthogonal to mom and o1, ie. o1 and o2  define plane.

    G4double x2, x1;

    do
    {
        x1=(G4UniformRand()-.5)*2;
        x2=(G4UniformRand()-.5)*2;
    } while (sqr(x1) +sqr(x2) > 1.);                      /* Loop checking, 31.08.2015, G.Folger */ // or random is badly broken.....

    return G4ThreeVector(r*(x1*o1.unit() + x2*o2.unit() - 1.5* mom.unit()));



    /*
     * // Get a point uniformly distributed on the surface of a sphere,
     * // with z < 0.
     *   G4double b = r*G4UniformRand();  // impact parameter
     *   G4double phi = G4UniformRand()*2*pi;
     *   G4double x = b*std::cos(phi);
     *   G4double y = b*std::sin(phi);
     *   G4double z = -std::sqrt(r*r-b*b);
     *   z *= 1.001; // Get position a little bit out of the sphere...
     *   point.setX(x);
     *   point.setY(y);
     *   point.setZ(z);
     */
}

//----------------------------------------------------------------------------

void G4BinaryCascade::ClearAndDestroy(G4KineticTrackVector * ktv)
//----------------------------------------------------------------------------
{
    std::vector<G4KineticTrack *>::iterator i;
    for(i = ktv->begin(); i != ktv->end(); ++i)
        delete (*i);
    ktv->clear();
}

//----------------------------------------------------------------------------
void G4BinaryCascade::ClearAndDestroy(G4ReactionProductVector * rpv)
//----------------------------------------------------------------------------
{
    std::vector<G4ReactionProduct *>::iterator i;
    for(i = rpv->begin(); i != rpv->end(); ++i)
        delete (*i);
    rpv->clear();
}

//----------------------------------------------------------------------------
void G4BinaryCascade::PrintKTVector(G4KineticTrackVector * ktv, std::string comment)
//----------------------------------------------------------------------------
{
    if (comment.size() > 0 ) G4cout << "G4BinaryCascade::PrintKTVector() " << comment << G4endl;
    if (ktv) {
        G4cout << "  vector: " << ktv << ", number of tracks: " << ktv->size()
			         << G4endl;
        std::vector<G4KineticTrack *>::iterator i;
        G4int count;

        for(count = 0, i = ktv->begin(); i != ktv->end(); ++i, ++count)
        {
            G4KineticTrack * kt = *i;
            G4cout << "  track n. " << count;
            PrintKTVector(kt);
        }
    } else {
        G4cout << "G4BinaryCascade::PrintKTVector():No KineticTrackVector given " << G4endl;
    }
}
//----------------------------------------------------------------------------
void G4BinaryCascade::PrintKTVector(G4KineticTrack * kt, std::string comment)
//----------------------------------------------------------------------------
{
    if (comment.size() > 0 ) G4cout << "G4BinaryCascade::PrintKTVector() "<< comment << G4endl;
    if ( kt ){
        G4cout << ", id: " << kt << G4endl;
        G4ThreeVector pos = kt->GetPosition();
        G4LorentzVector mom = kt->Get4Momentum();
        G4LorentzVector tmom = kt->GetTrackingMomentum();
        const G4ParticleDefinition * definition = kt->GetDefinition();
        G4cout << "    definition: " << definition->GetPDGEncoding() << " pos: "
                << 1/fermi*pos << " R: " << 1/fermi*pos.mag() << " 4mom: "
                << 1/MeV*mom <<"Tr_mom" <<  1/MeV*tmom << " P: " << 1/MeV*mom.vect().mag()
                << " M: " << 1/MeV*mom.mag() << G4endl;
        G4cout <<"    trackstatus: "<<kt->GetState() << " isParticipant " << (kt->IsParticipant()?"T":"F")   <<G4endl;
    } else {
        G4cout << "G4BinaryCascade::PrintKTVector(): No Kinetictrack given" << G4endl;
    }
}


//----------------------------------------------------------------------------
G4double G4BinaryCascade::GetIonMass(G4int Z, G4int A)
//----------------------------------------------------------------------------
{
    G4double mass(0);
    if ( Z > 0 && A >= Z )
    {
        mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z,A);

    } else if ( A > 0 && Z>0 )
    {
        // charge Z > A; will happen for light nuclei with pions involved.
        mass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(A,A);

    } else if ( A >= 0 && Z<=0 )
    {
        // all neutral, or empty nucleus
        mass = A * G4Neutron::Neutron()->GetPDGMass();

    } else if ( A == 0  )
    {
        // empty nucleus, except maybe pions
        mass = 0;

    } else
    {
        G4cerr << "G4BinaryCascade::GetIonMass() - invalid (A,Z) = ("
                << A << "," << Z << ")" <<G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade::GetIonMass() - giving up");

    }
	 //  G4cout << "G4BinaryCascade::GetIonMass() Z, A, mass " << Z << " " << A << " " << mass << G4endl;
    return mass;
}
G4ReactionProductVector * G4BinaryCascade::FillVoidNucleusProducts(G4ReactionProductVector * products)
{
    // return product when nucleus is destroyed, i.e. charge=0, or theTargetList.size()=0
    G4double Esecondaries(0.);
    G4LorentzVector psecondaries;
    std::vector<G4KineticTrack *>::iterator iter;
    std::vector<G4ReactionProduct *>::iterator rpiter;
    decayKTV.Decay(&theFinalState);

    for(iter = theFinalState.begin(); iter != theFinalState.end(); ++iter)
    {
        G4ReactionProduct * aNew = new G4ReactionProduct((*iter)->GetDefinition());
        aNew->SetMomentum((*iter)->Get4Momentum().vect());
        aNew->SetTotalEnergy((*iter)->Get4Momentum().e());
        aNew->SetCreatorModel(theBIC_ID);
        Esecondaries +=(*iter)->Get4Momentum().e();
        psecondaries +=(*iter)->Get4Momentum();
        aNew->SetNewlyAdded(true);
        //G4cout << " Particle Ekin " << aNew->GetKineticEnergy() << G4endl;
        products->push_back(aNew);
    }

    // pull out late particles from collisions
    //theCollisionMgr->Print();
    while(theCollisionMgr->Entries() > 0)        /* Loop checking, 31.08.2015, G.Folger */
    {
        G4CollisionInitialState *
        collision = theCollisionMgr->GetNextCollision();

        if ( ! collision->GetTargetCollection().size() ){
            G4KineticTrackVector * lates = collision->GetFinalState();
            if ( lates->size() == 1 ) {
                G4KineticTrack * atrack=*(lates->begin());
                //PrintKTVector(atrack, " late particle @ void Nucl ");

                G4ReactionProduct * aNew = new G4ReactionProduct(atrack->GetDefinition());
                aNew->SetMomentum(atrack->Get4Momentum().vect());
                aNew->SetTotalEnergy(atrack->Get4Momentum().e());
                // FIXME: should take creator model from atrack:
                // aNew->SetCreatorModel(atrack->GetCreatorModel());
                Esecondaries +=atrack->Get4Momentum().e();
                psecondaries +=atrack->Get4Momentum();
                aNew->SetNewlyAdded(true);
                products->push_back(aNew);

            }
        }
        theCollisionMgr->RemoveCollision(collision);

    }

    // decay must be after loop on Collisions, and Decay() will delete entries in theSecondaryList, refered
    //   to by Collisions.
    decayKTV.Decay(&theSecondaryList);

    // Correct for momentum transfered to Nucleus
    G4ThreeVector transferCorrection(0);
    if ( (theSecondaryList.size() + theCapturedList.size()) > 0)
    {
    	transferCorrection= theMomentumTransfer /(theSecondaryList.size() + theCapturedList.size());
    }

    for(iter = theSecondaryList.begin(); iter != theSecondaryList.end(); ++iter)
    {
        G4ReactionProduct * aNew = new G4ReactionProduct((*iter)->GetDefinition());
        (*iter)->Update4Momentum((*iter)->Get4Momentum().vect()+transferCorrection);
        aNew->SetMomentum((*iter)->Get4Momentum().vect());
        aNew->SetTotalEnergy((*iter)->Get4Momentum().e());
        aNew->SetCreatorModel(theBIC_ID);
        Esecondaries +=(*iter)->Get4Momentum().e();
        psecondaries +=(*iter)->Get4Momentum();
        if ( (*iter)->IsParticipant() ) aNew->SetNewlyAdded(true);
        products->push_back(aNew);
    }


    for(iter = theCapturedList.begin(); iter != theCapturedList.end(); ++iter)
    {
        G4ReactionProduct * aNew = new G4ReactionProduct((*iter)->GetDefinition());
        (*iter)->Update4Momentum((*iter)->Get4Momentum().vect()+transferCorrection);
        aNew->SetMomentum((*iter)->Get4Momentum().vect());
        aNew->SetTotalEnergy((*iter)->Get4Momentum().e());
        aNew->SetCreatorModel(theBIC_ID);
        Esecondaries +=(*iter)->Get4Momentum().e();
        psecondaries +=(*iter)->Get4Momentum();
        aNew->SetNewlyAdded(true);
        products->push_back(aNew);
    }

    G4double SumMassNucleons(0.);
    G4LorentzVector pNucleons(0.);
    for(iter = theTargetList.begin(); iter != theTargetList.end(); ++iter)
    {
        SumMassNucleons += (*iter)->GetDefinition()->GetPDGMass();
        pNucleons += (*iter)->Get4Momentum();
    }

    G4double Ekinetic=theProjectile4Momentum.e() + initial_nuclear_mass - Esecondaries - SumMassNucleons;
     #ifdef debug_BIC_FillVoidnucleus
        G4LorentzVector deltaP=theProjectile4Momentum + G4LorentzVector(initial_nuclear_mass) -
    		                psecondaries - pNucleons;
        //G4cout << "BIC::FillVoidNucleus() nucleons : "<<theTargetList.size() << " ,  T: " << Ekinetic <<
    	//	     ", deltaP " <<  deltaP << " deltaPNoNucl " << deltaP + pNucleons << G4endl;
     #endif
    if (Ekinetic > 0. && theTargetList.size()){
        Ekinetic /= theTargetList.size();
    } else {
        G4double Ekineticrdm(0);
        if (theTargetList.size()) Ekineticrdm = ( 0.1 + G4UniformRand()*5.) * MeV;	// leave some  Energy for Nucleons
      	G4double TotalEkin(Ekineticrdm);
        for (rpiter=products->begin(); rpiter!=products->end(); ++rpiter){
        	TotalEkin+=(*rpiter)->GetKineticEnergy();
        }
        G4double correction(1.);
        if ( std::abs(Ekinetic) < 20*perCent * TotalEkin ){
        	correction=1. + (Ekinetic-Ekineticrdm)/TotalEkin;   // Ekinetic < 0 == IS < FS, need to reduce energies
        }
        #ifdef debug_G4BinaryCascade
			else {
				G4cout << "BLIC::FillVoidNucleus() fail correction, Ekinetic, TotalEkin " << Ekinetic << ""<< TotalEkin << G4endl;
			}
        #endif

        for (rpiter=products->begin(); rpiter!=products->end(); ++rpiter){
	    	(*rpiter)->SetKineticEnergy((*rpiter)->GetKineticEnergy()*correction);  // this sets kinetic & total energy
        	(*rpiter)->SetMomentum((*rpiter)->GetTotalMomentum() * (*rpiter)->GetMomentum().unit());

        }

        Ekinetic=Ekineticrdm*correction;
        if (theTargetList.size())Ekinetic /= theTargetList.size();

	}

    for(iter = theTargetList.begin(); iter != theTargetList.end(); ++iter) {
    	// set Nucleon it to be hit - as it is in fact
    	(*iter)->Hit();
        G4ReactionProduct * aNew = new G4ReactionProduct((*iter)->GetDefinition());
        aNew->SetKineticEnergy(Ekinetic);
        aNew->SetMomentum(aNew->GetTotalMomentum() * ((*iter)->Get4Momentum().vect().unit()));
        aNew->SetNewlyAdded(true);
        aNew->SetCreatorModel(theBIC_ID);
        products->push_back(aNew);
        Esecondaries += aNew->GetTotalEnergy();
        psecondaries += G4LorentzVector(aNew->GetMomentum(),aNew->GetTotalEnergy() );
    }
    psecondaries=G4LorentzVector(0);
    for (rpiter=products->begin(); rpiter!=products->end(); ++rpiter){
    	psecondaries += G4LorentzVector((*rpiter)->GetMomentum(),(*rpiter)->GetTotalEnergy() );
    }



    G4LorentzVector initial4Mom=theProjectile4Momentum + G4LorentzVector(initial_nuclear_mass);

    //G4cout << "::FillVoidNucleus()final e/p conservation initial" <<initial4Mom
    //	<< " final " << psecondaries << " delta " << initial4Mom-psecondaries << G4endl;

    G4ThreeVector SumMom=psecondaries.vect();

    SumMom=initial4Mom.vect()-SumMom;
    G4int loopcount(0);

    std::vector<G4ReactionProduct *>::reverse_iterator reverse;  // start to correct last added first
    while ( SumMom.mag() > 0.1*MeV && loopcount++ < 10)           /* Loop checking, 31.08.2015, G.Folger */
    {
	 	G4int index=products->size();
		for (reverse=products->rbegin(); reverse!=products->rend(); ++reverse, --index){
			SumMom=initial4Mom.vect();
			for (rpiter=products->begin(); rpiter!=products->end(); ++rpiter){
				SumMom-=(*rpiter)->GetMomentum();
			}

			G4double p=((*reverse)->GetMomentum()).mag();
			(*reverse)->SetMomentum(  p*(((*reverse)->GetMomentum()+SumMom).unit()));

		}
    }


    return products;
}
G4ReactionProductVector * G4BinaryCascade::HighEnergyModelFSProducts(G4ReactionProductVector * products,
        G4KineticTrackVector * secondaries)
{
    std::vector<G4KineticTrack *>::iterator iter;
    for(iter = secondaries->begin(); iter != secondaries->end(); ++iter)
    {
        G4ReactionProduct * aNew = new G4ReactionProduct((*iter)->GetDefinition());
        aNew->SetMomentum((*iter)->Get4Momentum().vect());
        aNew->SetTotalEnergy((*iter)->Get4Momentum().e());
        aNew->SetNewlyAdded(true);
        // FixMe: should take creator model from atrack:
        // aNew->SetCreatorModel(atrack->GetCreatorModel());

        //G4cout << " Particle Ekin " << aNew->GetKineticEnergy() << G4endl;
        products->push_back(aNew);
    }
    const G4ParticleDefinition* fragment = 0;
    if (currentA == 1 && currentZ == 0) {
        fragment = G4Neutron::NeutronDefinition();
    } else if (currentA == 1 && currentZ == 1) {
        fragment = G4Proton::ProtonDefinition();
    } else if (currentA == 2 && currentZ == 1) {
        fragment = G4Deuteron::DeuteronDefinition();
    } else if (currentA == 3 && currentZ == 1) {
        fragment = G4Triton::TritonDefinition();
    } else if (currentA == 3 && currentZ == 2) {
        fragment = G4He3::He3Definition();
    } else if (currentA == 4 && currentZ == 2) {
        fragment = G4Alpha::AlphaDefinition();;
    } else {
        fragment =
                G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(currentZ,currentA,0.0);
    }
    if (fragment != 0) {
        G4ReactionProduct * theNew = new G4ReactionProduct(fragment);
        theNew->SetMomentum(G4ThreeVector(0,0,0));
        theNew->SetTotalEnergy(massInNucleus);
        // FixMe: should take creator model from ???:
        // aNew->SetCreatorModel(???->GetCreatorModel());
        //theNew->SetFormationTime(??0.??);
        //G4cout << " Nucleus (" << currentZ << ","<< currentA << "), mass "<< massInNucleus << G4endl;
        products->push_back(theNew);
    }
    return products;
}

void G4BinaryCascade::PrintWelcomeMessage()
{
    G4cout <<"Thank you for using G4BinaryCascade. "<<G4endl;
}

//----------------------------------------------------------------------------
void G4BinaryCascade::DebugApplyCollisionFail(G4CollisionInitialState * collision,
      G4KineticTrackVector * products)
{
   G4bool havePion=false;
   if (products)
   {
      for ( std::vector<G4KineticTrack *>::iterator i =products->begin(); i != products->end(); i++)
      {
         G4int PDGcode=std::abs((*i)->GetDefinition()->GetPDGEncoding());
         if (std::abs(PDGcode)==211 || PDGcode==111 ) havePion=true;
      }
   }
   if ( !products  || havePion)
   {
		const G4BCAction &action= *collision->GetGenerator();
      G4cout << " Collision " << collision << ", type: "<< typeid(action).name()
							                    << ", with NO products! " <<G4endl;
      G4cout << G4endl<<"Initial condition are these:"<<G4endl;
      G4cout << "proj: "<<collision->GetPrimary()->GetDefinition()->GetParticleName()<<G4endl;
      PrintKTVector(collision->GetPrimary());
      for(size_t it=0; it<collision->GetTargetCollection().size(); it++)
      {
         G4cout << "targ: "
               <<collision->GetTargetCollection()[it]->GetDefinition()->GetParticleName()<<G4endl;
      }
      PrintKTVector(&collision->GetTargetCollection(),std::string(" Target particles"));
   }
   //  if ( lateParticleCollision ) G4cout << " Added late particle--------------------------"<<G4endl;
   //  if ( lateParticleCollision && products ) PrintKTVector(products, " reaction products");
}

//----------------------------------------------------------------------------

G4bool G4BinaryCascade::CheckChargeAndBaryonNumber(G4String where)
{
   static G4int lastdA(0), lastdZ(0);
   G4int iStateA = the3DNucleus->GetMassNumber() + projectileA;
   G4int iStateZ = the3DNucleus->GetCharge()  + projectileZ;

   G4int fStateA(0);
   G4int fStateZ(0);

   std::vector<G4KineticTrack *>::iterator i;
   G4int CapturedA(0), CapturedZ(0);
   G4int secsA(0), secsZ(0);
   for ( i=theCapturedList.begin(); i!=theCapturedList.end(); ++i) {
      CapturedA += (*i)->GetDefinition()->GetBaryonNumber();
      CapturedZ += G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus);
   }

   for ( i=theSecondaryList.begin(); i!=theSecondaryList.end(); ++i) {
      if ( (*i)->GetState() != G4KineticTrack::inside ) {
         secsA += (*i)->GetDefinition()->GetBaryonNumber();
         secsZ += G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus);
      }
   }

   for ( i=theFinalState.begin(); i!=theFinalState.end(); ++i) {
      fStateA += (*i)->GetDefinition()->GetBaryonNumber();
      fStateZ += G4lrint((*i)->GetDefinition()->GetPDGCharge()/eplus);
   }

   G4int deltaA= iStateA -  secsA - fStateA -currentA - lateA;
   G4int deltaZ= iStateZ -  secsZ - fStateZ -currentZ - lateZ;

#ifdef debugCheckChargeAndBaryonNumberverbose	
	G4cout << where <<" A: iState= "<< iStateA<<", secs= "<< secsA<< ", fState= "<< fStateA<< ", current= "<<currentA<< ", late= " <<lateA << G4endl;   
	G4cout << where <<" Z: iState= "<< iStateZ<<", secs= "<< secsZ<< ", fState= "<< fStateZ<< ", current= "<<currentZ<< ", late= " <<lateZ << G4endl;   
#endif

   if (deltaA != 0  || deltaZ!=0 ) {
      if (deltaA != lastdA || deltaZ != lastdZ ) {
         G4cout << "baryon/charge imbalance - " << where << G4endl
               << "deltaA " <<deltaA<<", iStateA "<<iStateA<< ",  CapturedA "<<CapturedA <<",  secsA "<<secsA
               << ", fStateA "<<fStateA << ", currentA "<<currentA << ", lateA "<<lateA << G4endl
               << "deltaZ "<<deltaZ<<", iStateZ "<<iStateZ<< ",  CapturedZ "<<CapturedZ <<",  secsZ "<<secsZ
               << ", fStateZ "<<fStateZ << ", currentZ "<<currentZ << ", lateZ "<<lateZ << G4endl<< G4endl;
         lastdA=deltaA;
         lastdZ=deltaZ;
      }
   } else { lastdA=lastdZ=0;}

   return true;
}
//----------------------------------------------------------------------------
void G4BinaryCascade::DebugApplyCollision(G4CollisionInitialState * collision,
        G4KineticTrackVector * products)
{

    PrintKTVector(collision->GetPrimary(),std::string(" Primary particle"));
    PrintKTVector(&collision->GetTargetCollection(),std::string(" Target particles"));
    PrintKTVector(products,std::string(" Scatterer products"));

#ifdef dontUse
    G4double thisExcitation(0);
    //  excitation energy from this collision
    //  initial state:
    G4double initial(0);
    G4KineticTrack * kt=collision->GetPrimary();
    initial +=  kt->Get4Momentum().e();

    G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;

    initial +=  RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition());
    initial -=  RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding());
    G4cout << "prim. E/field/Barr/Sum " << kt->Get4Momentum().e()
        		          << " " << RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition())
        		          << " " << RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding())
        		          << " " << initial << G4endl;;

    G4KineticTrackVector ktv=collision->GetTargetCollection();
    for ( unsigned int it=0; it < ktv.size(); it++)
    {
        kt=ktv[it];
        initial +=  kt->Get4Momentum().e();
        thisExcitation += kt->GetDefinition()->GetPDGMass()
     		    		         - kt->Get4Momentum().e()
     		    		         - RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition());
        //     initial +=  RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition());
        //     initial -=  RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding());
        G4cout << "Targ. def/E/field/Barr/Sum " <<  kt->GetDefinition()->GetPDGEncoding()
  			          << " " << kt->Get4Momentum().e()
  			          << " " << RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition())
  			          << " " << RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding())
  			          << " " << initial <<" Excit " << thisExcitation << G4endl;;
    }

    G4double final(0);
    G4double mass_out(0);
    G4int product_barions(0);
    if ( products )
    {
        for ( unsigned int it=0; it < products->size(); it++)
        {
            kt=(*products)[it];
            final +=  kt->Get4Momentum().e();
            final +=  RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition());
            final +=  RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding());
            if ( kt->GetDefinition()->GetBaryonNumber()==1 ) product_barions++;
            mass_out += kt->GetDefinition()->GetPDGMass();
            G4cout << "sec. def/E/field/Barr/Sum " << kt->GetDefinition()->GetPDGEncoding()
  	    		         << " " << kt->Get4Momentum().e()
  	    		         << " " << RKprop->GetField(kt->GetDefinition()->GetPDGEncoding(),kt->GetPosition())
  	    		         << " " << RKprop->GetBarrier(kt->GetDefinition()->GetPDGEncoding())
  	    		         << " " << final << G4endl;;
        }
    }


    G4int finalA = currentA;
    G4int finalZ = currentZ;
    if ( products )
    {
        finalA -= product_barions;
        finalZ -= GetTotalCharge(*products);
    }
    G4double delta = GetIonMass(currentZ,currentA) - (GetIonMass(finalZ,finalA) + mass_out);
    G4cout << " current/final a,z " << currentA << " " << currentZ << " "<< finalA<< " "<< finalZ
            <<  " delta-mass " << delta<<G4endl;
    final+=delta;
    mass_out  = GetIonMass(finalZ,finalA);
    G4cout << " initE/ E_out/ Mfinal/ Excit " << currentInitialEnergy
            << " " <<   final << " "
            <<  mass_out<<" "
            <<  currentInitialEnergy - final - mass_out
            << G4endl;
    currentInitialEnergy-=final;
#endif
}

//----------------------------------------------------------------------------
G4bool G4BinaryCascade::DebugFinalEpConservation(const G4HadProjectile & aTrack,
        G4ReactionProductVector* products)
//----------------------------------------------------------------------------
{
    G4ReactionProductVector::iterator iter;
    G4double Efinal(0);
    G4ThreeVector pFinal(0);
    if (std::abs(theParticleChange.GetWeightChange() -1 ) > 1e-5 )
    {
        G4cout <<" BIC-weight change " << theParticleChange.GetWeightChange()<< G4endl;
    }

    for(iter = products->begin(); iter != products->end(); ++iter)
    {

           G4cout << " Secondary E - Ekin / p " <<
              (*iter)->GetDefinition()->GetParticleName() << " " <<
              (*iter)->GetTotalEnergy() << " - " <<
              (*iter)->GetKineticEnergy()<< " / " <<
              (*iter)->GetMomentum().x() << " " <<
              (*iter)->GetMomentum().y() << " " <<
              (*iter)->GetMomentum().z() << G4endl;
        Efinal += (*iter)->GetTotalEnergy();
        pFinal += (*iter)->GetMomentum();
    }

      G4cout << "e outgoing/ total : " << Efinal << " " << Efinal+GetFinal4Momentum().e()<< G4endl;
      G4cout << "BIC E/p delta " <<
            (aTrack.Get4Momentum().e()+theInitial4Mom.e() - Efinal)/MeV <<
            " MeV / mom " << (aTrack.Get4Momentum()  - pFinal ) /MeV << G4endl;

    return (aTrack.Get4Momentum().e() + theInitial4Mom.e() - Efinal)/aTrack.Get4Momentum().e() < perCent;

}
//----------------------------------------------------------------------------
G4bool G4BinaryCascade::DebugEpConservation(const G4String where)
//----------------------------------------------------------------------------
{
    G4cout << where << G4endl;
    G4LorentzVector psecs,    ptgts,    pcpts,    pfins;
    if (std::abs(theParticleChange.GetWeightChange() -1 ) > 1e-5 )
    {
        G4cout <<" BIC-weight change " << theParticleChange.GetWeightChange()<< G4endl;
    }

    std::vector<G4KineticTrack *>::iterator ktiter;
    for(ktiter = theSecondaryList.begin(); ktiter != theSecondaryList.end(); ++ktiter)
       {

              G4cout << " Secondary E - Ekin / p " <<
                 (*ktiter)->GetDefinition()->GetParticleName() << " " <<
                 (*ktiter)->Get4Momentum().e() << " - " <<
                 (*ktiter)->Get4Momentum().e() - (*ktiter)->Get4Momentum().mag() << " / " <<
                 (*ktiter)->Get4Momentum().vect() << G4endl;
           psecs += (*ktiter)->Get4Momentum();
       }

    for(ktiter = theTargetList.begin(); ktiter != theTargetList.end(); ++ktiter)
       {

              G4cout << " Target E - Ekin / p " <<
                 (*ktiter)->GetDefinition()->GetParticleName() << " " <<
                 (*ktiter)->Get4Momentum().e() << " - " <<
                 (*ktiter)->Get4Momentum().e() - (*ktiter)->Get4Momentum().mag() << " / " <<
                 (*ktiter)->Get4Momentum().vect() << G4endl;
           ptgts += (*ktiter)->Get4Momentum();
       }

    for(ktiter = theCapturedList.begin(); ktiter != theCapturedList.end(); ++ktiter)
        {

               G4cout << " Captured E - Ekin / p " <<
                  (*ktiter)->GetDefinition()->GetParticleName() << " " <<
                  (*ktiter)->Get4Momentum().e() << " - " <<
                  (*ktiter)->Get4Momentum().e() - (*ktiter)->Get4Momentum().mag() << " / " <<
                  (*ktiter)->Get4Momentum().vect() << G4endl;
            pcpts += (*ktiter)->Get4Momentum();
        }

    for(ktiter = theFinalState.begin(); ktiter != theFinalState.end(); ++ktiter)
        {

               G4cout << " Finals E - Ekin / p " <<
                  (*ktiter)->GetDefinition()->GetParticleName() << " " <<
                  (*ktiter)->Get4Momentum().e() << " - " <<
                  (*ktiter)->Get4Momentum().e() - (*ktiter)->Get4Momentum().mag() << " / " <<
                  (*ktiter)->Get4Momentum().vect() << G4endl;
            pfins += (*ktiter)->Get4Momentum();
        }

      G4cout << " Secondaries " << psecs << ", Targets " << ptgts << G4endl
    		  <<" Captured    " << pcpts << ", Finals  " << pfins << G4endl
    		  <<" Sum " << psecs + ptgts + pcpts + pfins << " PTransfer " << theMomentumTransfer
    		  <<" Sum+PTransfer " << psecs + ptgts + pcpts + pfins + theMomentumTransfer
    		  << G4endl<< G4endl;


    return true;

}

//----------------------------------------------------------------------------
G4ReactionProductVector * G4BinaryCascade::ProductsAddFakeGamma(G4ReactionProductVector *products )
//----------------------------------------------------------------------------
{
   //    else
//    {
//        G4ReactionProduct * aNew=0;
//        // return nucleus e and p
//        if  (fragment != 0 ) {
//            aNew = new G4ReactionProduct(G4Gamma::GammaDefinition());   // we only want to pass e/p
//            aNew->SetMomentum(fragment->GetMomentum().vect());
//            aNew->SetTotalEnergy(fragment->GetMomentum().e());
//            delete fragment;
//            fragment=0;
//        } else if (products->size() == 0) {
//            // FixMe GF: for testing without precompound, return 1 gamma of 0.01 MeV in +x
//#include "G4Gamma.hh"
//            aNew = new G4ReactionProduct(G4Gamma::GammaDefinition());
//            aNew->SetMomentum(G4ThreeVector(0.01*MeV,0,0));
//            aNew->SetTotalEnergy(0.01*MeV);
//        }
//        if ( aNew != 0 ) products->push_back(aNew);
//    }
   return products;
}

//----------------------------------------------------------------------------
