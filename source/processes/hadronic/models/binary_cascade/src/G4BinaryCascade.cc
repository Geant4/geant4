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
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4LorentzRotation.hh"
#include "G4BinaryCascade.hh"
#include "G4KineticTrackVector.hh"
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

#include "G4FermiPhaseSpaceDecay.hh"

#include <algorithm>
#include "G4ShortLivedConstructor.hh"
#include <typeinfo>

//   turn on general debugging info, and consistency checks
//#define debug_G4BinaryCascade 1

//  more detailed debugging -- deprecated  
//#define debug_1_BinaryCascade 1

//  specific debuuging info per method or functionality
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
//#define debug_BIC_FinalNucleusMomentum 1
//#define debug_BIC_FindFragments 1
//#define debug_BIC_BuildTargetList 1
//#define debug_BIC_FindCollision 1

//  C O N S T R U C T O R S   A N D   D E S T R U C T O R S
//

G4BinaryCascade::G4BinaryCascade() : 
G4VIntraNuclearTransportModel("Binary Cascade")
{
  // initialise the resonance sector
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  theCollisionMgr = new G4CollisionManager;
  theDecay=new G4BCDecay;
  theLateParticle= new G4BCLateParticle;
  theImR.push_back(theDecay);
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

  theExcitationHandler = new G4ExcitationHandler;
  SetDeExcitation(new G4PreCompoundModel(theExcitationHandler));
  SetMinEnergy(0.0*GeV);
  SetMaxEnergy(10.1*GeV);
  //PrintWelcomeMessage();
  thePrimaryEscape = true;
  thePrimaryType = 0;
}


G4BinaryCascade::G4BinaryCascade(const G4BinaryCascade& )
: G4VIntraNuclearTransportModel("Binary Cascade")
{
}


G4BinaryCascade::~G4BinaryCascade()
{
  ClearAndDestroy(&theTargetList);
  ClearAndDestroy(&theSecondaryList);
  ClearAndDestroy(&theCapturedList);
  ClearAndDestroy(&theProjectileList);
  delete thePropagator;
  delete theCollisionMgr;
  std::for_each(theImR.begin(), theImR.end(), Delete<G4BCAction>());
  delete theLateParticle;
  delete theExcitationHandler;
  delete theH1Scatterer;
}

//----------------------------------------------------------------------------

//
//      I M P L E M E N T A T I O N
//


//----------------------------------------------------------------------------
G4HadFinalState * G4BinaryCascade::ApplyYourself(const G4HadProjectile & aTrack,
//----------------------------------------------------------------------------
							G4Nucleus & aNucleus)
{
  static G4int eventcounter=0;
  
//   if ( eventcounter == 0 ) {
//      SetEpReportLevel(3);   // report non conservation with model etc.
//      G4double relativeLevel = 1*perCent;
//      G4double absoluteLevel = 2*MeV;
//      SetEnergyMomentumCheckLevels(relativeLevel,absoluteLevel); 
//   }
  
  //if(eventcounter == 100*(eventcounter/100) )
  eventcounter++;
  if(getenv("BCDEBUG") ) G4cerr << " ######### Binary Cascade Reaction number starts ######### "<<eventcounter<<G4endl;

  G4LorentzVector initial4Momentum = aTrack.Get4Momentum();
  G4ParticleDefinition * definition = const_cast<G4ParticleDefinition *>(aTrack.GetDefinition());

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
      G4cerr << "You are using G4BinaryCascade for projectiles other than nucleons or pions."<<G4endl;
      G4cerr << "If you want to continue, please switch on the developer environment: "<<G4endl;
      G4cerr << "setenv I_Am_G4BinaryCascade_Developer 1 "<<G4endl<<G4endl;
      throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCascade - used for unvalid particle type - Fatal");
    }
  }

// keep primary
  thePrimaryType = definition;
  thePrimaryEscape = false;

// try until an interaction will happen
  G4ReactionProductVector * products = 0;
  G4int interactionCounter = 0;
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

    the3DNucleus->Init(aNucleus.GetA_asInt(), aNucleus.GetZ_asInt());
    thePropagator->Init(the3DNucleus);
    //      GF Leak on kt??? but where to delete?
    G4KineticTrack * kt;// = new G4KineticTrack(definition, 0., initialPosition, initial4Momentum);
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
      products = Propagate(secondaries, the3DNucleus);
    } while(! products );  // until we FIND a collision...

    if(++interactionCounter>99) break;
  } while(products->size() == 0);  // ...untill we find an ALLOWED collision

  if(products->size()==0)
  {
    if(getenv("BCDEBUG") ) G4cerr << " ######### Binary Cascade Reaction number void ######### "<<eventcounter<<G4endl;
    theParticleChange.SetStatusChange(isAlive);
    theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
    theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    return &theParticleChange;
  }
//  G4cout << "BIC Applyyourself: number of products " << products->size() << G4endl;

// Fill the G4ParticleChange * with products
  theParticleChange.SetStatusChange(stopAndKill);
  G4ReactionProductVector::iterator iter;

  for(iter = products->begin(); iter != products->end(); ++iter)
  {
    G4DynamicParticle * aNew =
      new G4DynamicParticle((*iter)->GetDefinition(),
			    (*iter)->GetTotalEnergy(),
			    (*iter)->GetMomentum());
    // FixMe: should I use "position" or "time" specifyed AddSecondary() methods?
    theParticleChange.AddSecondary(aNew);

  }

  // DebugEpConservation(track, products);
              
  ClearAndDestroy(products);
  delete products;

  delete the3DNucleus;
  the3DNucleus = NULL;  // protect from wrong usage...

  if(getenv("BCDEBUG") ) G4cerr << " ######### Binary Cascade Reaction number ends ######### "<<eventcounter<<G4endl;
  	if (std::abs(theParticleChange.GetWeightChange() -1 ) > 1e-5 )
	{
	   G4cout <<" BIC-fin-weight change " << theParticleChange.GetWeightChange()<< G4endl;
	}
  return &theParticleChange;
}

//----------------------------------------------------------------------------
G4ReactionProductVector * G4BinaryCascade::Propagate(
//----------------------------------------------------------------------------
		G4KineticTrackVector * secondaries, G4V3DNucleus * nucleus)
{
  G4ping debug("debug_G4BinaryCascade");
  debug.push_back("trial");
  debug.dump();
#ifdef debug_BIC_Propagate
   G4cout << "G4BinaryCascade Propagate starting -------------------------------------------------------" <<G4endl;
#endif

   // *GF* FIXME ? in propagate mode this test is wrong! Could be in Apply....
  if(nucleus->GetMassNumber() == 1) // 1H1 is special case
  {
      #ifdef debug_BIC_Propagate
	  G4cout << " special case 1H1.... " << G4endl;
      #endif
     return Propagate1H1(secondaries,nucleus);
  }

  G4ReactionProductVector * products = new G4ReactionProductVector;
  the3DNucleus = nucleus;
  theOuterRadius = the3DNucleus->GetOuterRadius();
  theCurrentTime=0;
// build theSecondaryList, theProjectileList and theCapturedList
  ClearAndDestroy(&theCapturedList);
  ClearAndDestroy(&theSecondaryList);
  theSecondaryList.clear();
  ClearAndDestroy(&theProjectileList);
  ClearAndDestroy(&theFinalState);
  std::vector<G4KineticTrack *>::iterator iter;

  BuildTargetList();

   #ifdef debug_BIC_GetExcitationEnergy
     G4cout << "ExcitationEnergy0 " << GetExcitationEnergy() << G4endl;
   #endif

  thePropagator->Init(the3DNucleus);


  theCutOnP=90*MeV;
  if(nucleus->GetMass()>30) theCutOnP = 70*MeV;
  if(nucleus->GetMass()>60) theCutOnP = 50*MeV;
  if(nucleus->GetMass()>120) theCutOnP = 45*MeV;

  G4double StartingTime=DBL_MAX;        // Search for minimal formation time Uzhi
  for(iter = secondaries->begin(); iter != secondaries->end(); ++iter)    // Uzhi
  {                                                                       // Uzhi
    if((*iter)->GetFormationTime() < StartingTime)                        // Uzhi
        StartingTime = (*iter)->GetFormationTime();                       // Uzhi
  }                                                                       // Uzhi

  for(iter = secondaries->begin(); iter != secondaries->end(); ++iter)
  {
//  G4cout << " Formation time : " << (*iter)->GetDefinition()->GetParticleName() << " " 
// 	 << (*iter)->GetFormationTime() << G4endl;
    G4double FormTime = (*iter)->GetFormationTime() - StartingTime;       // Uzhi
    (*iter)->SetFormationTime(FormTime);                                  // Uzhi
    if( (*iter)->GetState() == G4KineticTrack::undefined ) 
    {
       FindLateParticleCollision(*iter);
    } else
    {
       theSecondaryList.push_back(*iter);
#ifdef debug_BIC_Propagate
       G4cout << " Adding initial secondary " << *iter 
                              << " time" << (*iter)->GetFormationTime()
			      << ", state " << (*iter)->GetState() << G4endl;
#endif
    }
//    theCollisionMgr->Print();
    theProjectileList.push_back(new G4KineticTrack(*(*iter)));
  }
  FindCollisions(&theSecondaryList);
  secondaries->clear(); // Don't leave "G4KineticTrack *"s in two vectors
  delete secondaries; 

// if called stand alone, build theTargetList and find first collisions

  if(theCollisionMgr->Entries() == 0 )      //late particles ALWAYS create Entries
  {
	//G4cout << " no collsions -> return 0" << G4endl;
      delete products;
      return 0;
  }

// end of initialization: do the job now
// loop untill there are no more collisions


  G4bool haveProducts = false;
  G4int collisionCount=0;
  theMomentumTransfer=G4ThreeVector(0,0,0);
  while(theCollisionMgr->Entries() > 0)
  {
    // G4cout << "Propagate Captured size: " << theCapturedList.size() << G4endl;
    // FixMe: at the moment I absorb pi with mom < theCutOnp, and I capture
    //        p and n with mom < theCutOnP. What about antip, k and sigmas
    //        with mom < theCutOnP?
    if(Absorb()) {  // absorb secondaries
      haveProducts = true;
      //G4cout << "Absorb sucess " << G4endl;
    }

    if(Capture()) { // capture secondaries (nucleons)
      haveProducts = true;
      //G4cout << "Capture sucess " << G4endl;
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
       debug.push_back("======>    test 1"); debug.dump();
       if (!DoTimeStep(nextCollision->GetCollisionTime()-theCurrentTime) )
       {
	   // Check if nextCollision is still valid, ie. particle did not leave nucleus
	   if (theCollisionMgr->GetNextCollision() != nextCollision )
	   {
	       nextCollision = 0;
	   }
	}
       debug.push_back("======>    test 2"); debug.dump();
//       G4cerr <<"post- DoTimeStep 1"<<G4endl;

	if( nextCollision )
	{
           debug.push_back("======>    test 3"); debug.dump();
	   if (ApplyCollision(nextCollision))
	   {
              //G4cerr << "ApplyCollision sucess " << G4endl;
 	      haveProducts = true;
	      collisionCount++;
              debug.push_back("======>    test 4.1"); debug.dump();
           } else {
              //G4cerr << "ApplyCollision failure " << G4endl;
	      theCollisionMgr->RemoveCollision(nextCollision);
              debug.push_back("======>    test 4.2"); debug.dump();
           }
	}
        debug.push_back("======>    test 5"); debug.dump();
//       G4cerr <<"post-post- DoTimeStep 1"<<G4endl;
    }
  }
  
//--------- end of while on Collsions  

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
	//G4cout << " no products return 0" << G4endl;
    return products;
  }


#ifdef debug_BIC_Propagate
   G4cout << " Momentum transfer to Nucleus " << theMomentumTransfer << " " << theMomentumTransfer.mag() << G4endl;
   G4cout << "  Stepping particles out...... " << G4endl;
#endif

  StepParticlesOut();
  
  
  if ( theSecondaryList.size() > 0 )
  {
#ifdef debug_G4BinaryCascade
      G4cerr << "G4BinaryCascade: Warning, have active particles at end" << G4endl;
#endif
//  add left secondaries to FinalSate
       std::vector<G4KineticTrack *>::iterator iter;
       for ( iter =theSecondaryList.begin(); iter != theSecondaryList.end(); ++iter)
       {
	   theFinalState.push_back(*iter);
       }
       theSecondaryList.clear();

  }
  while ( theCollisionMgr->Entries() > 0 )
  {
#ifdef debug_G4BinaryCascade
     G4cerr << " Warning: remove left over collision " << G4endl;
#endif
      theCollisionMgr->RemoveCollision(theCollisionMgr->GetNextCollision());
  }

#ifdef debug_BIC_Propagate_Excitation

  PrintKTVector(&theProjectileList,std::string(" theProjectileList"));
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
     } while ( ++ntry < maxtry && ExcitationEnergy < 0 );
  }

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
// 	if ( ExcitationEnergy < 0. )
 	{
//#ifdef debug_G4BinaryCascade
//  	  G4cerr << "G4BinaryCascade-Warning: negative excitation energy ";
//  	  G4cerr <<ExcitationEnergy<<G4endl;
// 	   PrintKTVector(&theFinalState,std::string("FinalState"));
// 	  PrintKTVector(&theCapturedList,std::string("captured"));
// 	  G4cout << "negative ExE:Final 4Momentum .mag: " << GetFinal4Momentum()
// 	          << " "<< GetFinal4Momentum().mag()<< G4endl
// 	          << "negative ExE:FinalNucleusMom  .mag: " << GetFinalNucleusMomentum()
// 		  << " "<< GetFinalNucleusMomentum().mag()<< G4endl;
//#endif
	}
	ClearAndDestroy(products);
		//G4cout << "  negative Excitation E return empty products " << products << G4endl;
	return products;   // return empty products
  }

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
       if(fragment->GetA() >1.5)                                          // Uzhi
       {
	 if (theDeExcitation)                // pre-compound
	 {
	      // G4cout << " going to preco with fragment 4 mom " << fragment->GetMomentum() << G4endl;
             pFragment=fragment->GetMomentum();
             precompoundProducts= theDeExcitation->DeExcite(*fragment);
             delete fragment;
             fragment=0;
	 } else if (theExcitationHandler)    // de-excitation
	 {
	      // G4cout << " going to de-excit with fragment 4 mom " << fragment->GetMomentum() << G4endl;
             pFragment=fragment->GetMomentum();
	     precompoundProducts=theExcitationHandler->BreakItUp(*fragment);
             delete fragment;
             fragment=0;
	 }
       } else 
       {                                   // fragment->GetA() < 1.5
	 precompoundProducts = new G4ReactionProductVector();
	 std::vector<G4KineticTrack *>::iterator i;
         if ( theTargetList.size() == 1 )
	 {
             i=theTargetList.begin();
	     G4ReactionProduct * aNew = new G4ReactionProduct((*i)->GetDefinition());
	     aNew->SetTotalEnergy((*i)->GetDefinition()->GetPDGMass());       
	     aNew->SetMomentum(G4ThreeVector(0));// see boost for preCompoundProducts below..
	     precompoundProducts->push_back(aNew);
	 } 

         if ( theCapturedList.size() == 1 )                               // Uzhi
	 {                                                                // Uzhi
             i=theCapturedList.begin();                                   // Uzhi
	     G4ReactionProduct * aNew = new G4ReactionProduct((*i)->GetDefinition()); // Uzhi
	     aNew->SetTotalEnergy((*i)->GetDefinition()->GetPDGMass());   // Uzhi
	     aNew->SetMomentum(G4ThreeVector(0));// see boost below..     // Uzhi
	     precompoundProducts->push_back(aNew);                        // Uzhi
	 }                                                                // Uzhi
       }                            // End of fragment->GetA() < 1.5
  } else                            // End of if(fragment)
  {                                 // No fragment, can be neutrons only  // Uzhi
     precompoundProducts = new G4ReactionProductVector();                      
     
     if ( (theTargetList.size()+theCapturedList.size()) > 0 )
     {                                                                      
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
           eCMS=sumMass + (2*MeV*masses.size());     
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
             precompoundProducts->push_back(aNew);

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
             precompoundProducts->push_back(aNew);                           // Uzhi
             delete *aMom;                                                   // Uzhi
	  }                                                                  // Uzhi
	}                                                                    // Uzhi

	if(momenta) delete momenta;
     }  
  }                   // End if(!fragment)


  {
// fill in products the outgoing particles
     G4double Ekinout=0;
     G4LorentzVector pSumBic(0);
     size_t i(0);
     for(i = 0; i< theFinalState.size(); i++)
     {
       G4KineticTrack * kt = theFinalState[i];
       if(kt->GetDefinition()->IsShortLived())
       {
         G4KineticTrackVector * dec = kt->Decay();
	 G4KineticTrackVector::const_iterator jter;
 	 for(jter=dec->begin(); jter != dec->end(); jter++)
	 {
	   theFinalState.push_back(*jter);
	 }
	 delete dec;
       }
       else
       {
	 G4ReactionProduct * aNew = new G4ReactionProduct(kt->GetDefinition());
	 aNew->SetMomentum(kt->Get4Momentum().vect());
	 aNew->SetTotalEnergy(kt->Get4Momentum().e());
	 Ekinout += aNew->GetKineticEnergy();
	 pSumBic += kt->Get4Momentum();
	 if(kt->IsParticipant()) 
	 {
	   aNew->SetNewlyAdded(true);
	 }
	 else
	 {
	   aNew->SetNewlyAdded(false);
	 }
	 //G4cout << " Particle Ekin " << aNew->GetKineticEnergy() << G4endl;
	 products->push_back(aNew);

	 #ifdef debug_BIC_Propagate_finals
	 if (! kt->GetDefinition()->GetPDGStable() )
	 {
             if (kt->GetDefinition()->IsShortLived())
	     {
		G4cout << "final shortlived : ";
	     } else
	     {
		G4cout << "final non stable : ";
	     }
	     G4cout <<kt->GetDefinition()->GetParticleName()<< G4endl;
	 }
	 #endif
       }

     }
     //G4cout << " Total Ekin " << Ekinout << G4endl;
  }
// add precompound products to products
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
	 G4cout << " pProduct be4 boost " <<pProduct << G4endl;
#endif
	 pProduct *= precompoundLorentzboost;
#ifdef debug_BIC_Propagate_finals
	 G4cout << " pProduct aft boost " <<pProduct << G4endl;
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
  else
  {
   G4ReactionProduct * aNew=0;
// return nucleus e and p
   if  (fragment != 0 ) {
       aNew = new G4ReactionProduct(G4Gamma::GammaDefinition());   // we only want to pass e/p
       aNew->SetMomentum(fragment->GetMomentum().vect());
       aNew->SetTotalEnergy(fragment->GetMomentum().e());
       delete fragment;
       fragment=0;
   } else if (products->size() == 0) {
   // FixMe GF: for testing without precompound, return 1 gamma of 0.01 MeV in +x
#include "G4Gamma.hh"
       aNew = new G4ReactionProduct(G4Gamma::GammaDefinition());
       aNew->SetMomentum(G4ThreeVector(0.01*MeV,0,0));
       aNew->SetTotalEnergy(0.01*MeV);
   }
   if ( aNew != 0 ) products->push_back(aNew);
  }

//  G4cerr <<"mon - end of event       "<<G4endl;
  thePrimaryEscape = true;
	//G4cout << "  return products " << products << G4endl;
  return products;
}


//----------------------------------------------------------------------------
G4double G4BinaryCascade::GetExcitationEnergy()
//----------------------------------------------------------------------------
{

  G4ping debug("debug_ExcitationEnergy");
// get A and Z for the residual nucleus
  #if defined(debug_G4BinaryCascade) || defined(debug_BIC_GetExcitationEnergy)
     G4int finalA = theTargetList.size()+theCapturedList.size();
     G4int finalZ = GetTotalCharge(theTargetList)+GetTotalCharge(theCapturedList);
     if ( (currentA - finalA) != 0 || (currentZ - finalZ) != 0 )
     {
	G4cerr << "G4BIC:GetExcitationEnergy(): Nucleon counting error current/final{A,Z} " 
               << currentA << " " << finalA << " "<< currentZ << " " << finalZ << G4endl;
     }
  
  #endif

  G4double excitationE(0);
  G4double nucleusMass(0);
  if(currentZ>.5)
  {
     nucleusMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA);
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
  debug.push_back("====> current A, Z");
  debug.push_back(currentZ);
  debug.push_back(currentA);
  debug.push_back(finalZ);
  debug.push_back(finalA);
  debug.push_back(nucleusMass);
  debug.push_back(GetFinalNucleusMomentum().mag());
  debug.dump();
//  PrintKTVector(&theTargetList, std::string(" current target list info"));
  PrintKTVector(&theCapturedList, std::string(" current captured list info"));
  #endif

  excitationE = GetFinalNucleusMomentum().mag() - nucleusMass;

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
      initialExc = theInitial4Mom.mag()-
           G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);
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
  G4ParticleDefinition * definition;
  G4ThreeVector pos;
  G4LorentzVector mom;
// if there are nucleon hit by higher energy models, then SUM(momenta) != 0
  theInitial4Mom = G4LorentzVector(0,0,0,0);
  currentA=0;
  currentZ=0;
//  G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;
  while((nucleon = the3DNucleus->GetNextNucleon()) != NULL)
  {
// check if nucleon is hit by higher energy model.
     if ( ! nucleon->AreYouHit() )
     {
	definition = nucleon->GetDefinition();
	pos = nucleon->GetPosition();
	mom = nucleon->GetMomentum();
 //    G4cout << "Nucleus " << pos.mag()/fermi << " " << mom.e() << G4endl;
	theInitial4Mom += mom;
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
     massInNucleus = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA);
  } else if (currentZ==0 && currentA>=1 )
  {
     massInNucleus = currentA * G4Neutron::Neutron()->GetPDGMass();
  } else
  {
     G4cerr << "G4BinaryCascade::BuildTargetList(): Fatal Error - invalid nucleus (A,Z)=("
		<< currentA << "," << currentZ << ")" << G4endl;
     throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde::BuildTargetList()");
  }
  currentInitialEnergy=	theInitial4Mom.e();
  G4KineticTrackVector::iterator i;
  for(i = theProjectileList.begin() ; i != theProjectileList.end(); ++i)
  {
    currentInitialEnergy+= (*i)->GetTrackingMomentum().e();
  }
	
#ifdef debug_BIC_BuildTargetList
     G4cout << "G4BinaryCascade::BuildTargetList():  nucleus (A,Z)=("
		<< currentA << "," << currentZ << ") mass: " << massInNucleus <<
		", theInitial4Mom " << theInitial4Mom << 
		", currentInitialEnergy " << currentInitialEnergy << G4endl;
#endif		

}


//----------------------------------------------------------------------------
void  G4BinaryCascade::FindCollisions(G4KineticTrackVector * secondaries)
//----------------------------------------------------------------------------
{
  for(std::vector<G4KineticTrack *>::iterator i = secondaries->begin();
     i != secondaries->end(); ++i)
  {
#ifdef debug_G4BinaryCascade
    if ( (*i)->GetTrackingMomentum().mag2() < -1.*eV )
    {
      G4cout << "G4BinaryCascade::FindCollisions(): negative m2:" << (*i)->GetTrackingMomentum().mag2() << G4endl;
    }
#endif
     
    for(std::vector<G4BCAction *>::iterator j = theImR.begin();
        j!=theImR.end(); j++)
    {
//      G4cout << "G4BinaryCascade::FindCollisions: at action " << *j << G4endl; 
      const std::vector<G4CollisionInitialState *> & aCandList
          = (*j)->GetCollisions(*i, theTargetList, theCurrentTime);
      for(size_t count=0; count<aCandList.size(); count++)
      {
        theCollisionMgr->AddCollision(aCandList[count]);
//	G4cout << "====================== New Collision ================="<<G4endl;
//	theCollisionMgr->Print();
      }
    }
  }
}

//----------------------------------------------------------------------------
void  G4BinaryCascade::FindDecayCollision(G4KineticTrack * secondary)
//----------------------------------------------------------------------------
{
    if ( secondary->GetTrackingMomentum().mag2() < -1.*eV )
    {
      G4cout << "G4BinaryCascade::FindDecayCollision(): negative m2:" << secondary->GetTrackingMomentum().mag2() << G4endl;
    } 
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
    G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;
    if ( secondary->GetTrackingMomentum().mag2() < -1.*eV )
    {
      G4cout << "G4BinaryCascade::FindDecayCollision(): negative m2:" << secondary->GetTrackingMomentum().mag2() << G4endl;
    } 

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

    //  for barions, correct for fermi energy
    G4int PDGcode=std::abs(secondary->GetDefinition()->GetPDGEncoding());
    if ( PDGcode > 1000 ) 
    {
       G4double initial_Efermi = RKprop->GetField(G4Neutron::Neutron()->GetPDGEncoding(),
                                          secondary->GetPosition()); 
       secondary->Update4Momentum(secondary->Get4Momentum().e() - initial_Efermi);
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
  G4ping debug("debug_ApplyCollision");
#ifdef debug_BIC_ApplyCollision
  G4cerr << "G4BinaryCascade::ApplyCollision start"<<G4endl;
  theCollisionMgr->Print();
#endif
  G4KineticTrack * primary = collision->GetPrimary();
  G4KineticTrackVector target_collection=collision->GetTargetCollection();
  G4bool haveTarget=target_collection.size()>0;
  if( haveTarget && (primary->GetState() != G4KineticTrack::inside) )
  {
#ifdef debug_G4BinaryCascade
     G4cout << "G4BinaryCasacde::ApplyCollision(): StateError " << primary << G4endl;
     G4KineticTrackVector debug;
     debug.push_back(primary);
     PrintKTVector(&debug,std::string("primay- ..."));
     PrintKTVector(&target_collection,std::string("... targets"));
//*GF*     throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde::ApplyCollision()");
#endif
     return false;
  }

  G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;
 
#ifdef debug_BIC_ApplyCollision
      G4cout << "ApplyCollisions : projte 4mom " << primary->GetTrackingMomentum()<< G4endl;
#endif

  G4int initialBaryon(0);
  G4int initialCharge(0);
  G4LorentzVector mom4Primary(0);
  
  if (primary->GetState() == G4KineticTrack::inside)
  {
     initialBaryon = primary->GetDefinition()->GetBaryonNumber();
     initialCharge = G4lrint(primary->GetDefinition()->GetPDGCharge());
  }

  G4int PDGcode=std::abs(primary->GetDefinition()->GetPDGEncoding());
  mom4Primary=primary->Get4Momentum();

// for primary resonances, subtract neutron ( = proton) field ( ie. add std::abs(field))
  G4double initial_Efermi(0);
  if (primary->GetState() == G4KineticTrack::inside ) {
     initial_Efermi=RKprop->GetField(primary->GetDefinition()->GetPDGEncoding(),primary->GetPosition());

     if ( PDGcode > 1000 && PDGcode != 2112 && PDGcode != 2212 )
     {
	initial_Efermi = RKprop->GetField(G4Neutron::Neutron()->GetPDGEncoding(),
                                                   primary->GetPosition());
	primary->Update4Momentum(mom4Primary.e() - initial_Efermi);
     }

     std::vector<G4KineticTrack *>::iterator titer;
     for ( titer=target_collection.begin() ; titer!=target_collection.end(); ++titer)
     {
	G4ParticleDefinition * aDef=(*titer)->GetDefinition();
	G4int aCode=aDef->GetPDGEncoding();
	G4ThreeVector aPos=(*titer)->GetPosition();
	initial_Efermi+= RKprop->GetField(aCode, aPos);
     }
  }
//****************************************

  G4KineticTrackVector * products=0;
  products = collision->GetFinalState();

  G4bool lateParticleCollision= (!haveTarget) && products && products->size() == 1;

  #ifdef debug_BIC_ApplyCollision
        G4bool havePion=false;
        for ( std::vector<G4KineticTrack *>::iterator i =products->begin(); i != products->end(); i++)
   	{
       		G4int PDGcode=std::abs((*i)->GetDefinition()->GetPDGEncoding());
		if (std::abs(PDGcode)==211 || PDGcode==111 ) havePion=true;
	}	
     if ( !products  || havePion)
      {
            G4cout << " Collision " << collision << ", type: "<< typeid(*collision->GetGenerator()).name()
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
   #endif
//****************************************  

  // reset primary to initial state
  primary->Set4Momentum(mom4Primary);


  G4int lateBaryon(0), lateCharge(0);

  if ( lateParticleCollision )
  {  // for late particles, reset charges
        //G4cout << "lateP, initial B C state " << initialBaryon << " " 
        //        << initialCharge<< " " << primary->GetState() << " "<< primary->GetDefinition()->GetParticleName()<< G4endl;
      lateBaryon = initialBaryon;
      lateCharge = initialCharge;
      initialBaryon=initialCharge=0;
  }
  
  initialBaryon += collision->GetTargetBaryonNumber();
  initialCharge+=G4lrint(collision->GetTargetCharge()); 
  if(!products || products->size()==0 || !CheckPauliPrinciple(products))
  {
   #ifdef debug_BIC_ApplyCollision
     if (products) G4cout << " ======Failed Pauli =====" << G4endl;
     G4cerr << "G4BinaryCascade::ApplyCollision blocked"<<G4endl;
   #endif
     if (products) ClearAndDestroy(products);
     if ( ! haveTarget ) FindDecayCollision(primary);  // for decay, sample new decay
     delete products;
     return false;
  }

  if (primary->GetState() == G4KineticTrack::inside ) {   // if the primary was outside, nothing to correct
     G4double final_Efermi(0);
     G4KineticTrackVector resonances;
     for ( std::vector<G4KineticTrack *>::iterator i =products->begin(); i != products->end(); i++)
     {
	 G4int PDGcode=std::abs((*i)->GetDefinition()->GetPDGEncoding());
         //       G4cout << " PDGcode, state " << PDGcode << " " << (*i)->GetState()<<G4endl;
	 final_Efermi+=RKprop->GetField(PDGcode,(*i)->GetPosition());
	 if ( PDGcode > 1000 && PDGcode != 2112 && PDGcode != 2212 )
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
               ClearAndDestroy(products);
               if (target_collection.size() == 0 ) FindDecayCollision(primary);  // for decay, sample new decay
	       delete products;
	       return false;
	    }
            //	  G4cout << " correct resonance from /to " << mom.e() << " / " << newEnergy<< G4endl;
	    G4ThreeVector mom3=std::sqrt(newEnergy2 - mass2) * mom.vect().unit();
	    (*res)->Set4Momentum(G4LorentzVector(mom3,newEnergy));
	}
     }
  }

#ifdef debug_BIC_ApplyCollision
  DebugApplyCollision(collision, products);
#endif

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
          finalCharge+=G4lrint((*i)->GetDefinition()->GetPDGCharge());
       }
    } else {
       G4double tin=0., tout=0.; 
       if (((G4RKPropagation*)thePropagator)->GetSphereIntersectionTimes((*i),tin,tout))
       {
          if ( tin > 0 ) 
	  {
	     (*i)->SetState(G4KineticTrack::outside);
	  }
	  else if ( tout > 0 ) 
	  {
	     (*i)->SetState(G4KineticTrack::inside);
             finalBaryon+=(*i)->GetDefinition()->GetBaryonNumber();
             finalCharge+=G4lrint((*i)->GetDefinition()->GetPDGCharge());	     
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
       
       //G4cout << " PDGcode, state " << (*i)->GetDefinition()->GetPDGEncoding() << " " << (*i)->GetState()<<G4endl;

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
  if (primary) 
  {
     oldSecondaries.push_back(primary);
     primary->Hit();
  }

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

  debug.push_back("=======> we have hit nucleons <=======");
  
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
// ------ debug
    G4int count = 0;
// ------ end debug
    while(!CheckPauliPrinciple(products))
    {
// ------ debug
      count++;
// ------ end debug
      ClearAndDestroy(products);
      if(!absorber.FindProducts(*kt))
	throw G4HadronicException(__FILE__, __LINE__, 
	  "G4BinaryCascade::Absorb(): Cannot absorb a particle.");
    }
// ------ debug
//    G4cerr << "Absorb CheckPauliPrinciple count= " <<  count << G4endl;
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
			 << " " << capturedEnergy/particlesBelowCut << " " << 0.2*theCutOnP << G4endl;
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
  G4ParticleDefinition * definition;

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
  while( theSecondaryList.size() > 0 )
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
       G4cout << " NextCollision  * , Time= " << nextCollision << " "
       		<<timeToCollision<< G4endl; 
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
    if ( counter > 100 && theCollisionMgr->Entries() == 0)   // no collision, and stepping a while....
    {
        #ifdef debug_1_BinaryCascade
        PrintKTVector(&theSecondaryList,std::string("stepping 100 steps"));
	#endif
	FindCollisions(&theSecondaryList);
	counter=0;
	++countreset;
    }
  }
//  G4cerr <<"Finished capture loop "<<G4endl;

       //G4cerr <<"pre- DoTimeStep 4"<<G4endl;
  DoTimeStep(DBL_MAX);
       //G4cerr <<"post- DoTimeStep 4"<<G4endl;


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

  G4double s = pCM.mag2();
//  G4double m10 = massInNucleus; //pNucleus.mag();
  G4double m10 = GetIonMass(currentZ,currentA);
  G4double m20 = pFinals.mag();
  if( s-(m10+m20)*(m10+m20) < 0 )
  {
       #ifdef debug_BIC_CorrectFinalPandE
	G4cout << "G4BinaryCascade::CorrectFinalPandE() : error! " << G4endl;

	G4cout << "not enough mass to correct: mass, A,Z, mass(nucl), mass(finals) " 
              << std::sqrt(-s+(m10+m20)*(m10+m20)) << " " 
	      << currentA << " " << currentZ << " "
	      << m10 << " " << m20 
	      << G4endl;
	G4cerr << " -CorrectFinalPandE 4" << G4endl;

	PrintKTVector(&theFinalState," mass problem");
       #endif
      return;
  }

  // Three momentum in cm system
  G4double pInCM = std::sqrt((s-(m10+m20)*(m10+m20))*(s-(m10-m20)*(m10-m20))/(4.*s));
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
  // G4cout << "Entering ... "<<oldTarget<<G4endl;
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
       PrintKTVector(kt_gone_out, std::string("append to final state.."));
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
	 secondaryCharge_in += G4lrint((*iter)->GetDefinition()->GetPDGCharge());
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
	     currentZ -= G4lrint((*iter)->GetDefinition()->GetPDGCharge());
	     currentA -= (*iter)->GetDefinition()->GetBaryonNumber();
	   
	}
	   
     }
#ifdef debug_BIC_CorrectBarionsOnBoundary
   G4cout << " CorrectBarionsOnBoundary, aft, A, Z, sec-Z,A,m,m_in_nucleus "
       << currentA << " " << currentZ << " "
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
	 secondaryCharge_out += G4lrint((*iter)->GetDefinition()->GetPDGCharge());
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

//  G4cout << "CorrectBarionsOnBoundary,secondaryCharge_out, secondaryBarions_out"
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

     if (secondaries_out>1) correction /= secondaries_out;
#ifdef debug_BIC_CorrectBarionsOnBoundary
       G4cout << "DoTimeStep,currentZ,currentA,"
	      << "secondaries_out,"
              <<"secondaryCharge_out,secondaryBarions_out,"
	      <<"energy correction,m_secondry,m_nucl_init,m_nucl_final "
 	 << " "<< currentZ << " "<< currentA <<" "
	 << secondaries_out << " " 
	 << secondaryCharge_out<<" "<<secondaryBarions_out<<" "
	 << correction << " "
 	 << secondaryMass_out << " "
 	 << mass_initial << " "
	 << mass_final << " "
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
	     G4RKPropagation * RKprop=(G4RKPropagation *)thePropagator;
	     (*iter)->SetState(G4KineticTrack::captured);
	     // Undo correction for Colomb Barrier
	     G4double barrier=RKprop->GetBarrier((*iter)->GetDefinition()->GetPDGEncoding());
	     (*iter)->UpdateTrackingMomentum((*iter)->GetTrackingMomentum().e() - barrier); 
	     if ( kt_fail == 0 ) kt_fail=new G4KineticTrackVector;
	     kt_fail->push_back(*iter);   
	     currentZ += G4lrint((*iter)->GetDefinition()->GetPDGCharge());
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
        << G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA)
        << " " << massInNucleus -G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA)
	<< G4endl;
#endif
  }
  
  return kt_fail;
}

					   
//----------------------------------------------------------------------------

G4Fragment * G4BinaryCascade::FindFragments()
//----------------------------------------------------------------------------
{

  G4int a = theTargetList.size()+theCapturedList.size();
#ifdef debug_BIC_FindFragments
  G4cout << "target, captured, secondary: "
         << theTargetList.size() << " " 
	 << theCapturedList.size()<< " "
	 << theSecondaryList.size()
	 << G4endl;
#endif
  G4int zTarget = 0;
  G4KineticTrackVector::iterator i;
  for(i = theTargetList.begin(); i != theTargetList.end(); ++i)
  {
      if((*i)->GetDefinition()->GetPDGCharge() == eplus)
      {
         zTarget++;
      }
  }

  G4int zCaptured = 0;
  G4LorentzVector CapturedMomentum=0;
  for(i = theCapturedList.begin(); i != theCapturedList.end(); ++i)
  {
      CapturedMomentum += (*i)->Get4Momentum();
      if((*i)->GetDefinition()->GetPDGCharge() == eplus)
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
 *          << G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z, a)
 * 	 << " / " << GetFinal4Momentum().mag()
 * 	 << " difference "
 * 	 <<  GetFinal4Momentum().mag() - G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z, a)
 * 	 << G4endl;
 */
//
//  if(getenv("BCDEBUG") ) G4cerr << "Fragment A, Z "<< a <<" "<< z<<G4endl;
  if ( z < 1 ) return 0;

  G4int holes = the3DNucleus->GetMassNumber() - theTargetList.size();
  G4int excitons = theCapturedList.size();
#ifdef debug_BIC_FindFragments
   G4cout << "Fragment: a= " << a
 	 << " z= " << z
 	 << " particles= " <<  excitons
 	 << " Charged= " << zCaptured
 	 << " holes= " << holes
 	 << " excitE= " <<GetExcitationEnergy()
 	 << " Final4Momentum= " << GetFinalNucleusMomentum()
 	 << " capturMomentum= " << CapturedMomentum
 	 << G4endl;
#endif

  G4Fragment * fragment = new G4Fragment(a,z,GetFinalNucleusMomentum());
  fragment->SetNumberOfHoles(holes);

//GF  fragment->SetNumberOfParticles(excitons-holes);
  fragment->SetNumberOfParticles(excitons);
  fragment->SetNumberOfCharged(zCaptured);
  G4ParticleDefinition * aIonDefinition =
       G4ParticleTable::GetParticleTable()->FindIon(a,z,0,z);
  fragment->SetParticleDefinition(aIonDefinition);

  return fragment;
}

//----------------------------------------------------------------------------

G4LorentzVector G4BinaryCascade::GetFinal4Momentum()
//----------------------------------------------------------------------------
{
// the initial 3-momentum will differ from 0, if nucleus created by string model.
  G4LorentzVector final4Momentum = theInitial4Mom;
  G4KineticTrackVector::iterator i;
  for(i = theProjectileList.begin() ; i != theProjectileList.end(); ++i)
  {
    final4Momentum += (*i)->GetTrackingMomentum();
    //G4cout << "Initial state: "<<(*i)->Get4Momentum()<<G4endl;
  }
  for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
  {
    final4Momentum -= (*i)->Get4Momentum();
  }

  if((final4Momentum.vect()/final4Momentum.e()).mag()>1.0 && currentA > 0)
  {
#  ifdef debug_BIC_Final4Momentum
     G4cerr << G4endl;
     G4cerr << "G4BinaryCascade::GetFinal4Momentum - Fatal"<<G4endl;
     G4KineticTrackVector::iterator i;
     G4cerr <<" GetFinal4Momentum: Initial nucleus "<<theInitial4Mom<<G4endl;
     for(i = theProjectileList.begin() ; i != theProjectileList.end(); ++i)
     {
       G4cerr << " Initial state (get4M), (trackingM): "
            <<(*i)->Get4Momentum()<<", " << (*i)->GetTrackingMomentum() <<","
	    <<(*i)->GetDefinition()->GetParticleName()<<G4endl;
     }
     for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
     {
       G4cerr <<" Final state: "<<(*i)->Get4Momentum()<<(*i)->GetDefinition()->GetParticleName()<<G4endl;
     }
     G4cerr<< " Final4Momentum = "<<final4Momentum <<" "<<final4Momentum.m()<<G4endl;
     G4cerr <<" current A, Z = "<< currentA<<", "<<currentZ<<G4endl;
     G4cerr << G4endl;
#  endif

     final4Momentum=G4LorentzVector(0);
  }
  return final4Momentum;
}

//----------------------------------------------------------------------------
G4LorentzVector G4BinaryCascade::GetFinalNucleusMomentum()
//----------------------------------------------------------------------------
{
// return momentum of nucleus for use with precompound model; also keep transformation to
//   apply to precompoud products.

  G4LorentzVector CapturedMomentum=0;
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
	boost=G4ThreeVector(0);
	NucleusMomentum=G4LorentzVector(0);
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
    G4ParticleDefinition * aHTarg = G4Proton::ProtonDefinition();
    G4double mass = aHTarg->GetPDGMass();
    if (nucleus->GetCharge() == 0) aHTarg = G4Neutron::NeutronDefinition();
    mass = aHTarg->GetPDGMass();
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
    while(!done && tryCount++ <200)
    {
      if(secs)
      {
       std::for_each(secs->begin(), secs->end(), DeleteKineticTrack());
       delete secs;
      }
      secs = theH1Scatterer->Scatter(*(*secondaries).front(), aTarget);
      for(size_t ss=0; secs && ss<secs->size(); ss++)
      {
//        G4cout << "1H1 " << (*secs)[ss]->GetDefinition()->GetParticleName()
//	       << ", shortlived? "<< (*secs)[ss]->GetDefinition()->IsShortLived()<< G4endl;
        if((*secs)[ss]->GetDefinition()->IsShortLived()) done = true;
      }
//    G4cout << G4endl;
    }
    size_t current(0);
    for(current=0; secs && current<secs->size(); current++)
    {
      if((*secs)[current]->GetDefinition()->IsShortLived())
      {
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
    //G4cout << "Through loop"<<G4endl;
    delete secs;
    for(iter = theFinalState.begin(); iter != theFinalState.end(); ++iter)
    {
      G4KineticTrack * kt = *iter;
      G4ReactionProduct * aNew = new G4ReactionProduct(kt->GetDefinition());
      aNew->SetMomentum(kt->Get4Momentum().vect());
      aNew->SetTotalEnergy(kt->Get4Momentum().e());
      products->push_back(aNew);
      #ifdef debug_1_BinaryCascade
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
    } while (sqr(x1) +sqr(x2) > 1.);

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
  if (comment.size() > 0 ) G4cout << comment << G4endl;
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
}
//----------------------------------------------------------------------------
void G4BinaryCascade::PrintKTVector(G4KineticTrack * kt, std::string comment)
//----------------------------------------------------------------------------
{
  if (comment.size() > 0 ) G4cout << comment << G4endl;
    G4cout << ", id: " << kt << G4endl;
    G4ThreeVector pos = kt->GetPosition();
    G4LorentzVector mom = kt->Get4Momentum();
    G4LorentzVector tmom = kt->GetTrackingMomentum();
    G4ParticleDefinition * definition = kt->GetDefinition();
    G4cout << "    definition: " << definition->GetPDGEncoding() << " pos: "
	   << 1/fermi*pos << " R: " << 1/fermi*pos.mag() << " 4mom: "
	   << 1/MeV*mom <<"Tr_mom" <<  1/MeV*tmom << " P: " << 1/MeV*mom.vect().mag() 
	   << " M: " << 1/MeV*mom.mag() << G4endl;
    G4cout <<"    trackstatus: "<<kt->GetState()<<G4endl;
}

//----------------------------------------------------------------------------
G4bool G4BinaryCascade::CheckDecay(G4KineticTrackVector * products)
//----------------------------------------------------------------------------
{
  G4int A = the3DNucleus->GetMassNumber();
  G4int Z = the3DNucleus->GetCharge();

  G4FermiMomentum fermiMom;
  fermiMom.Init(A, Z);

  const G4VNuclearDensity * density=the3DNucleus->GetNuclearDensity();

  G4KineticTrackVector::iterator i;
  G4ParticleDefinition * definition;

// ------ debug
  G4bool myflag = true;
// ------ end debug
  for(i = products->begin(); i != products->end(); ++i)
  {
    definition = (*i)->GetDefinition();
    if((definition->GetParticleName() != "delta++" )&&
       (definition->GetParticleName() != "delta+" )&&
       (definition->GetParticleName() != "delta0" )&&
       (definition->GetParticleName() != "delta-" ))
       continue;
       G4ThreeVector pos = (*i)->GetPosition();
       G4double d = density->GetDensity(pos);
       G4double pFermi= fermiMom.GetFermiMomentum(d);
       G4LorentzVector mom = (*i)->Get4Momentum();
       G4LorentzRotation boost(mom.boostVector()); 
       G4ThreeVector pion3(227*MeV * mom.vect().unit()); // 227 is decay product in rest frame
       G4LorentzVector pion(pion3, std::sqrt(sqr(140*MeV) +pion3.mag()));
     // G4cout << "pi rest " << pion << G4endl;
       pion = boost * pion;
     // G4cout << "pi lab  " << pion << G4endl;
// ------ debug
//     G4cout << "p:[" << (1/MeV)*pion.x() << " " << (1/MeV)*pion.y() << " "
//   	 << (1/MeV)*pion.z() << "] |p3|:"
//   	 << (1/MeV)*pion.vect().mag() << " E:" << (1/MeV)*pion.t() << " ] m: "
//   	 << (1/MeV)*pion.mag() << "  pos[" 
//   	 << (1/fermi)*pos.x() << " "<< (1/fermi)*pos.y() << " "
//   	 << (1/fermi)*pos.z() << "] |Dpos|: "
//   	 <<  " Pfermi:"
//   	 << (1/MeV)*pFermi << G4endl;   
// ------ end debug
       
     if( pion.vect().mag() < pFermi )
     {
// ------ debug
       myflag = false;
// ------ end debug
    }
  }
  return myflag;
//  return true;
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
      
   } else if ( A == 0 && std::abs(Z)<2 )
   {
      // empty nucleus, except maybe pions
      mass = 0;
      
   } else
   {
      G4cerr << "G4BinaryCascade::GetIonMass() - invalid (A,Z) = ("
              << A << "," << Z << ")" <<G4endl;
      G4Exception("G4BinaryCascade::GetIonMass() - giving up");
   }
   return mass;
}

void G4BinaryCascade::PrintWelcomeMessage()
{
  G4cout <<"Thank you for using G4BinaryCascade. "<<G4endl;
}

//----------------------------------------------------------------------------
void G4BinaryCascade::DebugApplyCollision(G4CollisionInitialState * collision, 
                                          G4KineticTrackVector * products)
{

  G4KineticTrackVector debug1;
  debug1.push_back(collision->GetPrimary());
  PrintKTVector(&debug1,std::string(" Primary particle"));
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
  G4double delta = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA) 
                   - (G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(finalZ,finalA) 
		       + mass_out); 
  G4cout << " current/final a,z " << currentA << " " << currentZ << " "<< finalA<< " "<< finalZ 
        <<  " delta-mass " << delta<<G4endl;
  final+=delta;
    mass_out  = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(finalZ,finalA);
  G4cout << " initE/ E_out/ Mfinal/ Excit " << currentInitialEnergy
         << " " <<   final << " "
	 <<  mass_out<<" " 
	 <<  currentInitialEnergy - final - mass_out
	 << G4endl;
   currentInitialEnergy-=final;	 
#endif
}

//----------------------------------------------------------------------------
void G4BinaryCascade::DebugEpConservation(const G4HadProjectile & aTrack,
					  G4ReactionProductVector* products)			   
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

//   G4cout << " Secondary E - Ekin / p " <<
//      (*iter)->GetDefinition()->GetParticleName() << " " <<
//      (*iter)->GetTotalEnergy() << " - " <<
//      (*iter)->GetKineticEnergy()<< " / " <<
//      (*iter)->GetMomentum().x() << " " <<
//      (*iter)->GetMomentum().y() << " " <<
//      (*iter)->GetMomentum().z() << G4endl;
      Efinal += (*iter)->GetTotalEnergy();
      pFinal += (*iter)->GetMomentum();
  }

//  G4cout << "e outgoing/ total : " << Efinal << " " << Efinal+GetFinal4Momentum().e()<< G4endl;
    G4cout << "BIC E/p delta " << 
    (aTrack.Get4Momentum().e()+ the3DNucleus->GetMass() - Efinal)/MeV <<
    " MeV / mom " << (aTrack.Get4Momentum()  - pFinal ) /MeV << G4endl;

}

//----------------------------------------------------------------------------
