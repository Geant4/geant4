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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
#include "G4ping.hh"

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

#include <algorithm>
#include "G4ShortLivedConstructor.hh"
#define debug_1_BinaryCascade 1
#define debug_G4BinaryCascade 1
//
//  C O N S T R U C T O R S   A N D   D E S T R U C T O R S
//

G4BinaryCascade::G4BinaryCascade() : G4VIntraNuclearTransportModel()
{
  // initialise the resonance sector
  G4ShortLivedConstructor ShortLived;
  ShortLived.ConstructParticle();

  theCollisionMgr = new G4CollisionManager;
  theScatterer = new G4Scatterer;
//  theScatterer = new G4ScattererStub;
  thePropagator = new G4RKPropagation;
  theCurrentTime = 0.;
  theCutOnP = 70*MeV;  //GF was 150..
  theCutOnPAbsorb= 0*MeV;
//  G4ExcitationHandler *
  theExcitationHandler = new G4ExcitationHandler;
  SetDeExcitation(new G4PreCompoundModel(theExcitationHandler));
  SetMinEnergy(0.0*GeV);
  SetMaxEnergy(10.1*GeV);
  PrintWelcomeMessage();
  thePrimaryEscape = true;
  thePrimaryType = 0;
}


G4BinaryCascade::G4BinaryCascade(const G4BinaryCascade& right)
{
}


G4BinaryCascade::~G4BinaryCascade()
{
  ClearAndDestroy(&theTargetList);
  ClearAndDestroy(&theSecondaryList);
  ClearAndDestroy(&theCapturedList);
  ClearAndDestroy(&theProjectileList);
  delete thePropagator;
  delete theScatterer;
  delete theCollisionMgr;
}

//----------------------------------------------------------------------------

//
//      I M P L E M E N T A T I O N
//


//----------------------------------------------------------------------------
G4VParticleChange * G4BinaryCascade::ApplyYourself(const G4Track & aTrack,
//----------------------------------------------------------------------------
							G4Nucleus & aNucleus)
{
  static G4int eventcounter=0;
  //if(eventcounter == 100*(eventcounter/100) )
  if(getenv("KCDEBUG") ) G4cerr << " ######### Reaction number starts ######### "<<eventcounter<<G4endl;
  eventcounter++;
  G4LorentzVector initial4Momentum = aTrack.GetDynamicParticle()->Get4Momentum();
  if(initial4Momentum.e()-initial4Momentum.m()<theCutOnP/2.)
  {
    return theDeExcitation->ApplyYourself(aTrack, aNucleus);
  }

  theParticleChange.Initialize(aTrack);
// initialize the G4V3DNucleus from G4Nucleus
  the3DNucleus = new G4Fancy3DNucleus;
  the3DNucleus->Init(aNucleus.GetN(), aNucleus.GetZ());

  thePropagator->Init(the3DNucleus);

// Build a KineticTrackVector with the G4Track
  G4KineticTrackVector * secondaries = new G4KineticTrackVector;
  G4ParticleDefinition * definition = aTrack.GetDefinition();
  G4ThreeVector initialPosition(0., 0., 0.); // will be set later

  if(!getenv("I_Am_G4BinaryCascade_Developer") )
  {
    if(definition!=G4Neutron::NeutronDefinition() &&
      definition!=G4Proton::ProtonDefinition() )
    {
      G4cerr << "You are using G4BinaryCascade for projectiles other than neutron and proton."<<G4endl;
      G4cerr << "If you want to continue, please switch on the developer environment: "<<G4endl;
      G4cerr << "setenv I_Am_G4BinaryCascade_Developer 1 "<<G4endl<<G4endl;
      G4Exception("G4BinaryCascade - used for unvalid particle type - Fatal");
    }
  }

// keep primary
  thePrimaryType = definition;
  thePrimaryEscape = false;

// try untill an interaction will happen
  G4ReactionProductVector * products = 0;
  G4double radius = the3DNucleus->GetOuterRadius()+3*fermi;
  do
  {
// reset status that could be changed in previous loop event
    theCollisionMgr->ClearAndDestroy();
    ClearAndDestroy(&theTargetList);  // clear and rebuild theTargetList
    BuildTargetList();
    G4KineticTrack * kt = new G4KineticTrack(definition, 0., initialPosition,
					     initial4Momentum);
// Note: secondaries has been cleared by Propagate() in the previous loop event
    secondaries->push_back(kt);

    while(theCollisionMgr->Entries() == 0)  // until we FIND a collision...
    {
      theCurrentTime=0;
      initialPosition=GetSpherePoint(1.1*radius, initial4Momentum);  // get random position
      kt->SetPosition(initialPosition);         // kt is in secondaries!!
      kt->SetState(G4KineticTrack::outside);
      FindCollisions(secondaries);
    }
    if(products != NULL)
    {  // free memory from previous loop event
      ClearAndDestroy(products);
      delete products;
    }
    products = Propagate(secondaries, the3DNucleus);

// --------- debug
// untill Propagate will be able to return some product...
//  if(1)
//    return &theParticleChange;
//  --------- end debug
  } while(products->size() == 0);  // ...untill we find an ALLOWED collision

//  G4cout << "HKM Applyyourself: number of products " << products->size() << G4endl;

// Fill the G4ParticleChange * with products
  G4int nProducts = products->size();
  theParticleChange.SetStatusChange(fStopAndKill);
  theParticleChange.SetNumberOfSecondaries(nProducts);
  G4ReactionProductVector::iterator iter;
  G4double Efinal=0;
  for(iter = products->begin(); iter != products->end(); ++iter)
  {
    G4DynamicParticle * aNew =
      new G4DynamicParticle((*iter)->GetDefinition(),
			    (*iter)->GetTotalEnergy(),
			    (*iter)->GetMomentum());
    if(getenv("KCDEBUG") ) 
    {
      if(abs(aNew->GetDefinition()->GetPDGEncoding()) >100
         && abs(aNew->GetDefinition()->GetPDGEncoding()) < 300) G4cout << "Pion info "<<aNew->GetDefinition()->GetPDGEncoding() <<" "<<aNew->GetKineticEnergy()<<G4endl;
    }
    // FixMe: should I use "position" or "time" specifyed AddSecondary() methods?
    theParticleChange.AddSecondary(aNew);

//   G4cout << " Secondary E - Ekin / p " <<
//      (*iter)->GetDefinition()->GetParticleName() << " " <<
//      (*iter)->GetTotalEnergy() << " - " <<
//      (*iter)->GetKineticEnergy()<< " / " <<
//      (*iter)->GetMomentum().x() << " " <<
//      (*iter)->GetMomentum().y() << " " <<
//      (*iter)->GetMomentum().z() << G4endl;
      Efinal += (*iter)->GetTotalEnergy();
  }

//  G4cout << "e outgoing/ total : " << Efinal << " " << Efinal+GetFinal4Momentum().e()<< G4endl;

  ClearAndDestroy(products);
  delete products;
  delete secondaries;

  delete the3DNucleus;
  the3DNucleus = NULL;  // protect from wrong usage...

  if(getenv("KCDEBUG") ) G4cerr << " ######### Reaction number ends ######### "<<eventcounter<<G4endl;
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
//  G4cerr <<"mon - Beginning of event       "<<G4endl;
   G4cout << "G4BinaryCascade Propagate starting -------------------------------------------------------" <<G4endl;
  G4ReactionProductVector * products = new G4ReactionProductVector;
  the3DNucleus = nucleus;
  theOuterRadius = the3DNucleus->GetOuterRadius();
  theCurrentTime=0;
// build theSecondaryList, theProjectileList and theCapturedList
  ClearAndDestroy(&theCapturedList);
  ClearAndDestroy(&theSecondaryList);
  ClearAndDestroy(&theProjectileList);
  ClearAndDestroy(&theFinalState);
  std::vector<G4KineticTrack *>::iterator iter;
  if(nucleus->GetMassNumber() == 1) // 1H1 is special case
  {
     return Propagate1H1(secondaries,nucleus);
  }

  for(iter = secondaries->begin(); iter != secondaries->end(); ++iter)
  {
    theSecondaryList.push_back(*iter);
    theProjectileList.push_back(new G4KineticTrack(*(*iter)));
  }
  secondaries->clear(); // Don't leave "G4KineticTrack *"s in two vectors

  thePropagator->Init(the3DNucleus);

// if called stand alone, build theTargetList and find first collisions

  if(theCollisionMgr->Entries() == 0)
  {
    the3DNucleus = nucleus;
    ClearAndDestroy(&theTargetList);
    BuildTargetList();
    FindCollisions(&theSecondaryList);
  }

//  thePropagator->Init(the3DNucleus);

// end of initialization: do the job now
// loop untill there are collisions


  G4bool haveProducts = false;
  G4int collisionCount=0;
  theMomentumTransfer=0;
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

       debug.push_back("======>    test 1"); debug.dump();
       if (!DoTimeStep(nextCollision->GetCollisionTime()-theCurrentTime) )
       {
	   // Check if nextCollision is still valid, ie. partcile did not leave nucleus
	   if (theCollisionMgr->GetNextCollision() != nextCollision )
	   {
	       nextCollision = 0;
	   }
	}
       debug.push_back("======>    test 2"); debug.dump();
//       G4cerr <<"post- DoTimeStep 1"<<G4endl;

	if( nextCollision )
	{
// GF for testing only, must be removed otherwise!
//	      theCollisionMgr->RemoveCollision(nextCollision);
           debug.push_back("======>    test 3"); debug.dump();
	   if (ApplyCollision(nextCollision))
// testing          if ( true ) 
	   {
// apply the collision
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

   debug.push_back("======>    #### through the first loop"); debug.dump();
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
    return products;
  }


#ifdef debug_G4BinaryCascade
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

// ------ debug

/*
 *   G4cout << "theProjectileList" << G4endl;
 *   PrintKTVector(&theProjectileList);
 *   G4cout << "theSecondaryList" << G4endl;
 *   PrintKTVector(&theSecondaryList);
 *   G4cout << "theTargetList size: " << theTargetList.size() << G4endl;
 * //  PrintKTVector(&theTargetList);
 *   G4cout << "theCapturedList" << G4endl;
 *   PrintKTVector(&theCapturedList);
 */

// ------ end debug

#ifdef debug_G4BinaryCascade
  G4cout << " ExcitE be4 Correct : " <<GetExcitationEnergy() << G4endl;
  G4cout << " Mom Transfered to nucleus : " << theMomentumTransfer << " " << theMomentumTransfer.mag() << G4endl;
#endif
PrintKTVector(&theFinalState,std::string(" FinalState uncorrected"));

  CorrectFinalPandE();
PrintKTVector(&theFinalState,std::string(" FinalState corrected"));

#ifdef debug_G4BinaryCascade
  G4cout << " ExcitE aft Correct : " <<GetExcitationEnergy() << G4endl;
#endif
//  CorrectFinalPandE();
//#ifdef debug_G4BinaryCascade
//  G4cout << " ExcitE aft2 Correct : " <<GetExcitationEnergy() << G4endl;
//#endif

//  G4cerr <<"mon - all pushed to limit 1"<<G4endl;
   G4double ExcitationEnergy=GetExcitationEnergy();
//  G4cerr <<"mon - all pushed to limit 2"<<G4endl;

//#ifdef HKM_DEBUG
//   G4cout << " Excitation Energy final, Ekinout, #collisions: "
//   << ExcitationEnergy << " "
//   << Ekinout << " "
//   << collisionCount <<G4endl;
//   G4cout << " Out from casc: " << theFinalState.size() << G4endl;
//#endif


  if ( ExcitationEnergy < 0. )
  {
// 	if ( ExcitationEnergy < 0. )
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
	}
	ClearAndDestroy(products);
	return products;   // return empty products
  }

//#ifdef HKM_DEBUG
//   G4cout << " Excitation Energy final, Ekinout, #collisions: "
//   << ExcitationEnergy << " "
//   << Ekinout << " "
//   << collisionCount <<G4endl;
//   G4cout << " Out from casc: " << theFinalState.size() << G4endl;
//#endif

// find a fragment and call the precompound model.
  G4Fragment * fragment = 0;
  G4ReactionProductVector * precompoundProducts = 0;
//  G4cerr <<"mon - entering deexcitat "<<G4endl;
   if ( ExcitationEnergy >= 0 ) // FixMe: GF temporary should we better re-start?
   {
//       G4Fragment *
       fragment = FindFragments();

    //  theDeExcitation =0;
       if(fragment && fragment->GetA() >1.5) // hpw @@@@ still to add the nucleon, if one exists.
       {
	 if (theDeExcitation)
	 {
             precompoundProducts= theDeExcitation->DeExcite(*fragment);
             delete fragment;
             fragment=0;
	 } else if (theExcitationHandler)
	 {
	     precompoundProducts=theExcitationHandler->BreakItUp(*fragment);
             delete fragment;
             fragment=0;
	 }
       } else {
         if ( theTargetList.size() == 1 )
	 {
	     precompoundProducts = new G4ReactionProductVector();
	     std::vector<G4KineticTrack *>::iterator i=theTargetList.begin();
	     G4ReactionProduct * aNew = new G4ReactionProduct((*i)->GetDefinition());
//	     G4LorentzVector final= GetFinalNucleusMomentum();
	     aNew->SetTotalEnergy((*i)->GetDefinition()->GetPDGMass());       //final.e());
//	     G4double mom(0);
//	     if ( final.e() > (*i)->GetDefinition()->GetPDGMass() )
//	     {
//	        mom=sqrt ( sqr(final.e()) - sqr((*i)->GetDefinition()->GetPDGMass()));
//	     }
	     aNew->SetMomentum(G4ThreeVector(0));		//mom*final.vect().unit());
	     precompoundProducts->push_back(aNew);
	  }
       }

   } 
#ifdef debug_G4BinaryCascade
   else {
       G4cerr << "Binary Cascade Error: negative Excitation Energy "<< G4endl;
   }
#endif
  {
// fill in products the outgoing particles
     G4double Ekinout=0;
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
         G4cout << " create product " << aNew << " from kt " << kt << G4endl;
	 aNew->SetMomentum(kt->Get4Momentum().vect());
	 aNew->SetTotalEnergy(kt->Get4Momentum().e());
	 Ekinout += aNew->GetKineticEnergy();
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
       }

     }
     //G4cout << " Total Ekin " << Ekinout << G4endl;
  }
// add precompound products to products
  if ( precompoundProducts )
  {
       std::vector<G4ReactionProduct *>::iterator j;
       for(j = precompoundProducts->begin(); j != precompoundProducts->end(); ++j)
       {
// boost back to system of moving nucleus
         G4LorentzVector pProduct((*j)->GetMomentum(),(*j)->GetTotalEnergy());
#ifdef debug_G4BinaryCascade
	 G4cout << " pProduct be4 boost " <<pProduct << G4endl;
#endif
	 pProduct *= precompoundLorentzboost;
#ifdef debug_G4BinaryCascade
	 G4cout << " pProduct aft boost " <<pProduct << G4endl;
#endif
         (*j)->SetTotalEnergy(pProduct.e());
         (*j)->SetMomentum(pProduct.vect());
	 (*j)->SetNewlyAdded(true);
	 products->push_back(*j);
       }

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
  return products;
}


//----------------------------------------------------------------------------
G4double G4BinaryCascade::GetExcitationEnergy()
//----------------------------------------------------------------------------
{

  G4ping debug("debug_ExcitationEnergy");
// get A and Z for the residual nucleus
  G4int finalA = theTargetList.size()+theCapturedList.size();
  G4int finalZ = GetTotalCharge(theTargetList)+GetTotalCharge(theCapturedList);

  G4double excitationE(0);
  G4double nucleusMass(0);
  if(currentZ>.5)
  {
     //G4cerr <<"mon - we know z           "<<currentZ<<" "<<currentA<<G4endl;
     nucleusMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA);
  } 
  else if (currentZ==0 && currentA==1 )
  {
     nucleusMass = G4Neutron::Neutron()->GetPDGMass();
  } 
  else
  {
     #ifdef debug_1_BinaryCascade
     G4cout << "G4BinaryCascade::GetExcitationEnergy(): Warning - invalid nucleus (A,Z)=("
	    << currentA << "," << currentZ << ")" << G4endl;
     #endif
     return 0;
  }

  debug.push_back("====> current A, Z");
  debug.push_back(currentZ);
  debug.push_back(currentA);
  debug.push_back(finalZ);
  debug.push_back(finalA);
  debug.push_back(nucleusMass);
  debug.push_back(GetFinalNucleusMomentum().mag());
  debug.dump();
  PrintKTVector(&theTargetList, std::string(" current target list info"));
  PrintKTVector(&theCapturedList, std::string(" current captured list info"));

  excitationE = GetFinalNucleusMomentum().mag() - nucleusMass;

#ifdef debug_G4BinaryCascade
// ------ debug
  if ( excitationE < 0 )
  {
     G4cout << "negative ExE final Ion mass" <<nucleusMass<< G4endl;
    if(finalZ>.5) G4cout << " Ecitation Energy, Finalnuclmom, nucl mass, excitE "
	               << GetFinalNucleusMomentum() << G4endl
	               <<excitationE << " "
	               << G4endl;

    if(finalZ>.5) G4cout << " final Excit : a,z, 4mom "
		         << finalA << " " << finalZ << " "
		<< GetFinalNucleusMomentum() <<G4endl;

    G4int A = the3DNucleus->GetMassNumber();
    G4int Z = the3DNucleus->GetCharge();
    G4double initialExc(0);
    if(Z>.5)
    {
      //G4cerr << "occation 1 "<< Z <<" "<<A<<G4endl;
      initialExc = theInitial4Mom.mag()-
           G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(Z, A);
    }

     //G4cout << " theInitial4Mom; 4.mag() " << theInitial4Mom << " "<< theInitial4Mom.mag() << G4endl;
     if(finalZ>.5)  G4cout << " finalNucleusMomentum; 4.mag() " << GetFinalNucleusMomentum() << " "<< GetFinalNucleusMomentum().mag() << G4endl;
  }

//   G4cout << "theCapturedList" << theCapturedList.size() << G4endl
//          << "theSecondaryList" << theSecondaryList.size() << G4endl;

// ------ end debug
#endif
  //  return excitationE > 0 ? excitationE : 0.0;
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
//   In the kinetic Model, the potential inside the nucleus is taken into account, and nucleons
//    are on mass shell.
	mom.setE( sqrt( mom.vect().mag2() + sqr(definition->GetPDGMass()) ) );
// 	if ( mom.e() + RKprop->GetField(definition->GetPDGEncoding(),pos) >definition->GetPDGMass() )
// 	{
// 	   G4cout << "G4BinaryCascade::BuildTargetList found free nucleon " <<G4endl
// 	   	  << "   4mom = " << mom << G4endl
// 		  << "    pos = " << pos << G4endl
// 		  << " mass, Field " <<definition->GetPDGMass() << " " <<RKprop->GetField(definition->GetPDGEncoding(),pos)
// 		  << G4endl;
// 	}


	G4KineticTrack * kt = new G4KineticTrack(definition, 0., pos, mom);
	kt->SetState(G4KineticTrack::inside);
	kt->SetNucleon(nucleon);
	theTargetList.push_back(kt);
	++currentA;
	if (definition->GetPDGCharge() > .5 ) ++currentZ;
     }
  }
  massInNucleus = 0;
  if(currentZ>.5)
  {
      //G4cerr << "occation 2 "<< currentZ <<" "<<currentA<<G4endl;
     massInNucleus = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(currentZ,currentA);
  } else if (currentZ==0 && currentA==1 )
  {
     massInNucleus = G4Neutron::Neutron()->GetPDGMass();
  } else
  {
     //G4cout << "G4BinaryCascade::BuildTargetList(): Warning - invalid nucleus (A,Z)=("
	//	<< currentA << "," << currentZ << ")" << G4endl;
  }

}


//----------------------------------------------------------------------------
void  G4BinaryCascade::FindCollisions(G4KineticTrackVector * secondaries)
//----------------------------------------------------------------------------
{
  G4KineticTrack * pkt;
  G4KineticTrack * tkt;
  G4double collisionTime;
/*
 * G4cout << " FindCollisions start" << endl;
 * theCollisionMgr->Print();
 */
  for(std::vector<G4KineticTrack *>::iterator i = secondaries->begin();
      i != secondaries->end(); ++i)
  {
    pkt = *i;
    // look for collisions with target particles
//      G4cerr << "G4BinaryCascade::ApplyCollision pre-collision time"
//             <<pkt->GetTrackingMomentum()<<" "<<pkt->GetTrackingMomentum().boostVector().mag()<<" " <<pkt->GetDefinition()->GetParticleName()
//	     <<G4endl;
    for(std::vector<G4KineticTrack *>::iterator j = theTargetList.begin();
	j != theTargetList.end(); ++j)
    {
      tkt = *j;
      collisionTime = theScatterer->GetTimeToInteraction(*pkt, *tkt);
      if(collisionTime == DBL_MAX)  // no collision
	continue;
      theCollisionMgr->AddCollision(collisionTime+theCurrentTime, pkt, tkt);
      G4cerr <<" !!!!!! debug collisions "<<collisionTime<<" "<<pkt->GetDefinition()->GetParticleName()<<G4endl;      
   }
    //  G4cerr << "G4BinaryCascade::ApplyCollision post-collision time"<<G4endl;
    // look for decay
    if(pkt->GetDefinition()->IsShortLived())
    {
      collisionTime = pkt->SampleResidualLifetime();
//      G4cerr <<" !!!!!! debug decay "<<collisionTime<<" "<<pkt->GetDefinition()->GetParticleName()<<G4endl;
      theCollisionMgr->AddCollision(collisionTime+theCurrentTime, pkt, 0);
/*
 *     G4cout << " shortlived: " << pkt->GetDefinition()->GetParticleName();
 *     G4cout << " shortlived: " << pkt->GetDefinition()->GetParticleName()
 *            << " coll.time (ns)=" << collisionTime/ns << G4endl;
 */
    }
  }

  G4cout << " FindCollisions found collisions" << endl;
  theCollisionMgr->Print();
 
}


//----------------------------------------------------------------------------
G4bool G4BinaryCascade::ApplyCollision(G4CollisionInitialState * collision)
//----------------------------------------------------------------------------
{
  G4ping debug("debug_ApplyCollision");
  //G4cerr << "G4BinaryCascade::ApplyCollision start"<<G4endl;
  G4KineticTrack * primary = collision->GetPrimary();
  G4KineticTrack * target = collision->GetTarget();


  // The projectile MUST not be outside of nucleus!
  if(primary->GetState() != G4KineticTrack::inside)
  {
     G4cout << "G4BinaryCasacde::ApplyCollision(): StateError " << primary << G4endl;
     G4KineticTrackVector debug;
     debug.push_back(primary);
     debug.push_back(target);
     PrintKTVector(&debug,std::string("primay- target"));
#ifdef debug_BinaryCascade
     G4Exception("G4BinaryCasacde::ApplyCollision()");
#else
     return false;
#endif
  }

  G4KineticTrackVector * products=0;

      G4cout << "ApplyCollisions : projte 4mom " << primary->GetTrackingMomentum()<< G4endl;
   if (target != 0 )
   {
      G4cout << "ApplyCollisions : target 4mom " << target->Get4Momentum()<< G4endl;
   }


  G4int initialBaryon = primary->GetDefinition()->GetBaryonNumber();
  G4int initialCharge(0);
  if(initialBaryon!=0) initialCharge+=static_cast<G4int>(primary->GetDefinition()->GetPDGCharge()+.1);
  if(target == NULL) // if target == 0 it is a decay
  {
       //G4cerr << "G4BinaryCascade::ApplyCollision pre-decay"<<G4endl;
    products = primary->Decay();
       //G4cerr << "G4BinaryCascade::ApplyCollision post-decay"<<G4endl;
     //G4cerr << " decaying a " << primary->GetDefinition()->GetParticleName() << G4endl;
     //G4cerr << "number of products created: " << products->size() << G4endl;
  }
  else
  {
    G4KineticTrack target_reloc(*target);
    initialBaryon += target->GetDefinition()->GetBaryonNumber();
    if(target->GetDefinition()->GetBaryonNumber()!=0) initialCharge+=static_cast<G4int>(target->GetDefinition()->GetPDGCharge()+.1);
  G4cerr << "G4BinaryCascade::ApplyCollision start 1"<<G4endl;
    products = theScatterer->Scatter(*primary, target_reloc);
  G4cerr << "G4BinaryCascade::ApplyCollision start 2"<<G4endl;

// @@@@@ also decay needs to be blocked...for after this release - hpw.

    if(!products || !CheckPauliPrinciple(products))
    {
       if (products) G4cout << " ======Failed Pauli =====" << G4endl;
       if (products) ClearAndDestroy(products);
       delete products;
       G4cerr << "G4BinaryCascade::ApplyCollision blocked"<<G4endl;
       return false;
    }
 }

  if(products == 0)
  {
     delete products;
     return false;
  }

// debug block
  #ifdef debug_1_BinaryCascade
  PrintKTVector(products,std::string(" Scatterer products"));
  #endif

//  G4cout << " ======Survive Pauli =====" << G4endl;
  G4int finalBaryon(0);
  G4int finalCharge(0);
  for(size_t ig=0;ig<products->size();ig++)
  {
    finalBaryon+=products->operator[](ig)->GetDefinition()->GetBaryonNumber();
    if(products->operator[](ig)->GetDefinition()->GetBaryonNumber()!=0) finalCharge+=static_cast<G4int>(products->operator[](ig)->GetDefinition()->GetPDGCharge()+.1);
// secondaries are created inside nucleus
    products->operator[](ig)->SetState(G4KineticTrack::inside);
  }

  currentA += finalBaryon-initialBaryon;
  currentZ += finalCharge-initialCharge;
  
  G4KineticTrackVector oldSecondaries;
  oldSecondaries.push_back(primary);
  G4cout << "Balancing debug ====== "
      <<finalBaryon<<" "<<initialBaryon<<" "<<finalCharge<<" "<<initialCharge<<G4endl;
  G4KineticTrackVector oldTarget;

  primary->Hit();
  if(target)
  {
    oldTarget.push_back(target);
    target->Hit();
  }

  debug.push_back("=======> we have hit nucleons <=======");
  
  UpdateTracksAndCollisions(&oldSecondaries, &oldTarget, products);
  ClearAndDestroy(&oldSecondaries);  // free memory of desappeared tracks
  ClearAndDestroy(&oldTarget);
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
      G4Exception("G4BinaryCascade::Absorb(): Cannot absorb a particle.");

    if(!absorber.FindProducts(*kt))
      G4Exception("G4BinaryCascade::Absorb(): Cannot absorb a particle.");

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
	G4Exception(
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
	 if( energy < theCutOnP )
	 {
	    capturedEnergy+=energy;
	    ++particlesBelowCut;
	 } else
	 {
	    ++particlesAboveCut;
	 }
     }
    }
  }
  if(particlesAboveCut==0 && particlesBelowCut>0 && capturedEnergy/particlesBelowCut<0.2*theCutOnP)
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
       G4double eFermi = sqrt( sqr(fermiMom.GetFermiMomentum(d)) + (*i)->Get4Momentum().mag2() );
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
  #ifdef debug_G4BinaryCascade
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
  G4double steplength=1.0 * fermi;
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
	if((kt->GetDefinition() == G4Proton::Proton()) ||                     // @@@ GF why only for nucleons?
	   (kt->GetDefinition() == G4Neutron::Neutron()))
        {
	  nsec++;
	  G4double tStep = steplength / ( kt->Get4Momentum().beta() * c_light );
	  G4cout << " minTimeStep, tStep Particle " <<minTimeStep << " " <<tStep
	         << " " <<kt->GetDefinition()->GetParticleName() << " 4mom " << kt->GetTrackingMomentum()<<G4endl;
	  if(tStep<minTimeStep)
	  {
	    minTimeStep = tStep;
//            G4cerr <<"Position "<<kt->GetPosition().mag()<<" "
//	           <<kt->GetTrackingMomentum().e()-kt->GetTrackingMomentum().mag()<<G4endl;
	  }
        }
      } else if ( kt->GetState() != G4KineticTrack::outside ){
          PrintKTVector(&theSecondaryList, std::string(" state ERROR....."));
          G4Exception("G4BinaryCascade::StepParticlesOut() particle not in nucleus");
      }
    }
//    G4cerr << "CaptureCount = "<<counter<<" "<<nsec<<" "<<minTimeStep<<" "<<1*ns<<G4endl;
    G4double timeToCollision=DBL_MAX;
    G4CollisionInitialState * nextCollision=0;
    if(theCollisionMgr->Entries() > 0)
    {
       nextCollision = theCollisionMgr->GetNextCollision();
       timeToCollision = nextCollision->GetCollisionTime()-theCurrentTime;
    }
    //G4cout << " minTimeStep, timeToCollisiont  " <<minTimeStep << " " <<timeToCollision << G4endl;
    if ( timeToCollision > minTimeStep )
    {
       //G4cerr <<"pre- DoTimeStep 2"<<G4endl;
	DoTimeStep(minTimeStep);
        ++counter;
       //G4cerr <<"post- DoTimeStep 2"<<G4endl;
    } else
    {
       // G4cerr <<"pre- DoTimeStep 3"<<G4endl;
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

    if(Capture(true))
    {
//       haveProducts = true;
       G4cout << "Capture sucess " << G4endl;
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
{
  if ( theFinalState.size() == 0 ) return;

  G4KineticTrackVector::iterator i;
  G4LorentzVector pNucleus=GetFinal4Momentum();
  G4LorentzVector pFinals(0);
  G4int nFinals(0);
  for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
  {
    pFinals += (*i)->Get4Momentum();
    ++nFinals;
#ifdef debug_G4BinaryCascade
    G4cout <<"CorrectFinalPandE a final " << (*i)->GetDefinition()->GetParticleName()
           << " 4mom " << (*i)->Get4Momentum()<< G4endl;
#endif
  }
#ifdef debug_G4BinaryCascade
  G4cout << "CorrectFinalPandE pN pF: " <<pNucleus << " " <<pFinals << G4endl;
#endif
  G4LorentzVector pCM=pNucleus + pFinals;

  G4LorentzRotation toCMS(-pCM.boostVector());
  pFinals *=toCMS;

#ifdef debug_G4BinaryCascade
  G4cout << "CorrectFinalPandE pCM, CMS pCM " << pCM << " " <<toCMS*pCM<< G4endl;
  G4cout << "CorrectFinal CMS pN pF " <<toCMS*pNucleus << " "
         <<pFinals << G4endl
	 <<" massInNucleus m(nucleus) m(finals) sqrt(s): " << massInNucleus << " " <<pNucleus.mag()<< " "
	 << pFinals.mag() << " " << pCM.mag() << G4endl;
#endif

  G4LorentzRotation toLab = toCMS.inverse();

  G4double s = pCM.mag2();
  G4double m10 = massInNucleus; //pNucleus.mag();
  G4double m20 = pFinals.mag();
  if( s-(m10+m20)*(m10+m20) < 0 )
  {
#ifdef debug_G4BinaryCascade
     G4cout << "G4BinaryCascade::CorrectFinalPandE() : error! " << G4endl;
#endif
      return;
  }

  // Three momentum in cm system
  G4double pInCM = sqrt((s-(m10+m20)*(m10+m20))*(s-(m10-m20)*(m10-m20))/(4.*s));
#ifdef debug_G4BinaryCascade
  G4cout <<" CorrectFinalPandE pInCM current/new : " <<(pFinals).vect().mag() << " " <<pInCM << G4endl;
#endif
  if ( pFinals.vect().mag() > pInCM )
  {
    G4ThreeVector p3finals=pInCM*pFinals.vect().unit();

//    G4ThreeVector deltap=(p3finals - pFinals.vect() ) / nFinals;
    G4double factor=pInCM/pFinals.vect().mag();
    G4LorentzVector qFinals(0);
    for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
    {
//      G4ThreeVector p3((toCMS*(*i)->Get4Momentum()).vect() + deltap);
      G4ThreeVector p3(factor*(toCMS*(*i)->Get4Momentum()).vect());
      G4LorentzVector p(p3,sqrt((*i)->Get4Momentum().mag2() + p3.mag2()));
      qFinals += p;
      p *= toLab;
#ifdef debug_G4BinaryCascade
      G4cout << " final p corrected: " << p << G4endl;
#endif
      (*i)->Set4Momentum(p);
    }
#ifdef debug_G4BinaryCascade
    G4cout << "CorrectFinalPandE nucleus corrected mass : " << GetFinal4Momentum() << " "
    		<<GetFinal4Momentum().mag() << G4endl
		<< " CMS pFinals , mag, 3.mag : " << qFinals << " " << qFinals.mag() << " " << qFinals.vect().mag()<< G4endl;
#endif
   }

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
  if(oldSecondaries != NULL)
  {
    if(!oldSecondaries->empty())
    {
      for(iter1 = oldSecondaries->begin(); iter1 != oldSecondaries->end();
	  ++iter1)
      {
	iter2 = std::find(theSecondaryList.begin(), theSecondaryList.end(),
			    *iter1);
	theSecondaryList.erase(iter2);
      }
      theCollisionMgr->RemoveTracksCollisions(oldSecondaries);
    }
  }

// remove old target from the target list
  if(oldTarget != NULL)
  {
    if(!oldTarget->empty())
    {
      for(iter1 = oldTarget->begin(); iter1 != oldTarget->end(); ++iter1)
      {
	iter2 = std::find(theTargetList.begin(), theTargetList.end(),
			    *iter1);
	theTargetList.erase(iter2);
      }
      theCollisionMgr->RemoveTracksCollisions(oldTarget);
    }
  }

  if(newSecondaries != NULL)
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
  G4ping debug("debug_G4BinaryCascade");
  debug.push_back("======> DoTimeStep 1"); debug.dump();
  G4bool success=true;
  G4cerr <<"G4BinaryCascade::DoTimeStep: enter "<< theTimeStep << G4endl;
  PrintKTVector(&theSecondaryList, std::string("DoTimeStep - theSecondaryList"));
//   PrintKTVector(&theTargetList, std::string("DoTimeStep - theTargetList"));
  std::vector<G4KineticTrack *>::iterator iter;
// Count particles in nucleus
  G4int secondaryBarions=0;
  G4int secondaryCharge=0;
  G4double secondaryMass=0;


  G4KineticTrackVector * kt_outside = new G4KineticTrackVector;
  for_each( theSecondaryList.begin(),theSecondaryList.end(),
           SelectFromKTV(kt_outside,G4KineticTrack::outside));
//  PrintKTVector(kt_outside, std::string("DoTimeStep - found outside"));	  

  G4KineticTrackVector * kt_inside = new G4KineticTrackVector;
  for_each( theSecondaryList.begin(),theSecondaryList.end(),
           SelectFromKTV(kt_inside, G4KineticTrack::inside));
//  PrintKTVector(kt_inside, std::string("DoTimeStep - found inside"));	  
//-----
    G4KineticTrackVector dummy;   // needed for re-usability
    G4cout << "NOW WE ARE ENTERING THE TRANSPORT"<<G4endl;

// =================== Here we move the particles  ===================  

     thePropagator->Transport(theSecondaryList, dummy, theTimeStep);
     
// =================== Here we move the particles  ===================  

//------
  PrintKTVector(&theSecondaryList, std::string("DoTimeStep - secondairies"));

     theMomentumTransfer += thePropagator->GetMomentumTransfer();
    debug.push_back("======> DoTimeStep 1.1.1"); debug.dump();
// Partclies which went INTO nucleus

  G4KineticTrackVector * kt_gone_in = new G4KineticTrackVector;
  for_each( kt_outside->begin(),kt_outside->end(),
           SelectFromKTV(kt_gone_in,G4KineticTrack::inside));
//  PrintKTVector(kt_gone_in, std::string("DoTimeStep - gone in"));	  

   //PrintKTVector(&theSecondaryList,std::string("aft trsprt....."));

// Partclies which  went OUT OF nucleus
  G4KineticTrackVector * kt_gone_out = new G4KineticTrackVector;
  for_each( kt_inside->begin(),kt_inside->end(),
           SelectFromKTV(kt_gone_out, G4KineticTrack::gone_out));

//  PrintKTVector(kt_gone_out, std::string("DoTimeStep - gone out"));	  

  CorrectBarionsOnBoundary(kt_gone_in,kt_gone_out); 

// Add track missing nucleus to addFinals
  for_each( kt_outside->begin(),kt_outside->end(),
           SelectFromKTV(kt_gone_out,G4KineticTrack::miss_nucleus));
  for_each( kt_outside->begin(),kt_outside->end(),
           SelectFromKTV(kt_gone_out,G4KineticTrack::gone_out));
    
PrintKTVector(kt_gone_out, std::string("pushing to final state.."));

// Partclies which could not leave nucleus,  captured...
  G4KineticTrackVector * kt_captured = new G4KineticTrackVector;
    for_each( theSecondaryList.begin(),theSecondaryList.end(),
           SelectFromKTV(kt_captured, G4KineticTrack::captured));

// Check no track is part in next collision, ie.
//  this step was to far, and collisions should not occur any more 

  if ( theCollisionMgr->Entries()> 0 )
  {
     if (kt_gone_out->size() )
     {
	G4KineticTrack * nextPrimary = theCollisionMgr->GetNextCollision()->GetPrimary();
	iter = find(kt_gone_out->begin(),kt_gone_out->end(),nextPrimary);
	if ( iter !=  kt_gone_out->end() )
	{
	   success=false;
	   G4cout << " DoTimeStep - WARNING: deleting current collision!" << G4endl;
	}
     }	
     if ( kt_captured->size() )
     {
	G4KineticTrack * nextPrimary = theCollisionMgr->GetNextCollision()->GetPrimary();
	iter = find(kt_captured->begin(),kt_captured->end(),nextPrimary);
	if ( iter !=  kt_captured->end() )
	{
	   success=false;
	   G4cout << " DoTimeStep - WARNING: deleting current collision!" << G4endl;
	}
     }	

  }
  UpdateTracksAndCollisions(kt_gone_out,0 ,0);


  if ( kt_captured->size() )
  {
     theCapturedList.insert(theCapturedList.end(),
                            kt_captured->begin(),kt_captured->end());
     UpdateTracksAndCollisions(kt_captured, NULL, NULL);
  }
  
#ifdef debug_G4BinaryCascade
  delete kt_inside;
  kt_inside = new G4KineticTrackVector;
  for_each( theSecondaryList.begin(),theSecondaryList.end(),
           SelectFromKTV(kt_inside, G4KineticTrack::inside));
   if ( currentZ != (GetTotalCharge(theTargetList) 
                    + GetTotalCharge(theCapturedList)
		    + GetTotalCharge(*kt_inside)) )
   {
      G4cout << " error-DoTimeStep, aft, A, Z, sec-Z,A,m, "
       << currentA << " "
       << currentZ << " "
       << GetTotalCharge(theTargetList)
          + GetTotalCharge(theCapturedList)
          + GetTotalCharge(*kt_inside) << " "
       << secondaryCharge << " "
       << secondaryBarions << " "
       << secondaryMass << " "
       << GetTotalCharge(theTargetList) << " " 
       << GetTotalCharge(theCapturedList) << " "
       << GetTotalCharge(*kt_inside) << " "
       << G4endl;
   }    
#endif

  delete kt_inside;
  delete kt_outside;
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

  if (in->size())
  {
     G4int secondaryBarions_in(0);
     G4int secondaryCharge_in(0);
     G4double secondaryMass_in(0);

     for ( iter =in->begin(); iter != in->end(); ++iter)
     {
	 if((*iter)->GetDefinition()->GetBaryonNumber()==1) 
	 {
	   secondaryCharge_in += static_cast<G4int>((*iter)->GetDefinition()->GetPDGCharge() + 0.1);
	   ++secondaryBarions_in;
	   secondaryMass_in += (*iter)->GetDefinition()->GetPDGMass();
	 }
     }
     G4double mass_initial= G4ParticleTable::GetParticleTable()->GetIonTable()->
		      GetIonMass(currentZ,currentA);
		      
     currentZ += secondaryCharge_in;
     currentA += secondaryBarions_in;
     
     G4double mass_final= G4ParticleTable::GetParticleTable()->GetIonTable()->
		      GetIonMass(currentZ,currentA);
     
     G4double correction= secondaryMass_in + mass_initial - mass_final;
     if (secondaryBarions_in>1) 
       {correction /= secondaryBarions_in;}

#ifdef debug_G4BinaryCascade
       G4cout << "DoTimeStep,currentZ,currentA,"
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
	if((*iter)->GetDefinition()->GetBaryonNumber()==1) 
	{
	   (*iter)->Update4Momentum((*iter)->GetTrackingMomentum().e() + correction);
	}
     }
#ifdef debug_G4BinaryCascade
   G4cout << " DoTimeStep, aft, A, Z, sec-Z,A,m,m_in_nucleus "
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
     G4int secondaryBarions_out(0);
     G4int secondaryCharge_out(0);
     G4double secondaryMass_out(0);

     for ( iter =out->begin(); iter != out->end(); ++iter)
     {
	 if((*iter)->GetDefinition()->GetBaryonNumber()==1) 
	 {
	   secondaryCharge_out += static_cast<G4int>((*iter)->GetDefinition()->GetPDGCharge() + 0.1);
	   ++secondaryBarions_out;
	   secondaryMass_out += (*iter)->GetDefinition()->GetPDGMass();
	 }
     }

     G4double mass_initial= G4ParticleTable::GetParticleTable()->GetIonTable()->
		              GetIonMass(currentZ,currentA);
     currentA -=secondaryBarions_out;
     currentZ -=secondaryCharge_out;

     if (currentA < 0 || currentZ < 0 || currentA < currentZ) 
       { G4Exception("G4BinaryCascade::CorrectBarionsOnBoundary() - fatal error");}

     G4double mass_final(0);
     G4double correction(0);

     if(currentZ>=1) 
     {
        mass_final = G4ParticleTable::GetParticleTable()->GetIonTable()->
		       GetIonMass(currentZ,currentA);
     } else if (currentA >=1 ){
        mass_final= currentA * G4Neutron::Neutron()->GetPDGMass();
     }	// else no nucleus left -> mass=0

     correction = mass_initial - mass_final - secondaryMass_out;
     if (secondaryBarions_out>1) correction /= secondaryBarions_out;
#ifdef debug_G4BinaryCascade
       G4cout << "DoTimeStep,currentZ,currentA,"
              <<"secondaryCharge_out,secondaryBarions_out,"
	      <<"energy correction,m_secondry,m_nucl_init,m_nucl_final "
 	 << " "<< currentZ << " "<< currentA <<" "
	 << secondaryCharge_out<<" "<<secondaryBarions_out<<" "
	 << correction << " "
 	 << secondaryMass_out << " "
 	 << mass_initial << " "
	 << mass_final << " "
	 << G4endl;
#endif
     PrintKTVector(out,std::string("out be4 correction"));
 
     for ( iter = out->begin(); iter != out->end(); ++iter)
     {
	if((*iter)->GetDefinition()->GetBaryonNumber()==1) 
	{
	   (*iter)->Update4Momentum((*iter)->GetTrackingMomentum().e() + correction);
	}
	theFinalState.push_back(*iter);
     }

     PrintKTVector(out,std::string("out AFTER correction"));

#ifdef debug_G4BinaryCascade
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
  G4cout << "Findfragment exitation E: " << GetExcitationEnergy() << G4endl;


  G4int a = theTargetList.size()+theCapturedList.size();
  G4cout << "target, captured, secondary: "
         << theTargetList.size() << " " 
	 << theCapturedList.size()<< " "
	 << theSecondaryList.size()
	 << G4endl;
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
//  G4cout <<"Findfragments- Captured size: " <<theCapturedList.size()<< G4endl;
  for(i = theCapturedList.begin(); i != theCapturedList.end(); ++i)
  {
//   G4cout << "Findfragments- Captured particle mom: " <<(*i)->Get4Momentum() << G4endl;
      CapturedMomentum += (*i)->Get4Momentum();
      if((*i)->GetDefinition()->GetPDGCharge() == eplus)
      {
	 zCaptured++;
      }
  }

  G4int z = zTarget+zCaptured;

  if ( z != (GetTotalCharge(theTargetList) + GetTotalCharge(theCapturedList)) )
  {
      G4cout << " FindFragment Counting error z a " << z << " " <<a << " "  
      << GetTotalCharge(theTargetList) << " " <<  GetTotalCharge(theCapturedList)<<
      G4endl;
      PrintKTVector(&theTargetList, std::string("theTargetList"));
      PrintKTVector(&theCapturedList, std::string("theCapturedList"));
  }
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
//  if(getenv("KCDEBUG") ) G4cerr << "Fragment A, Z "<< a <<" "<< z<<G4endl;
  if ( z < 1 ) return 0;
  G4Fragment * fragment = new G4Fragment(a,z,GetFinalNucleusMomentum());

  G4int holes = the3DNucleus->GetMassNumber() - theTargetList.size();
  fragment->SetNumberOfHoles(holes);

  G4int excitons = theCapturedList.size();
//GF  fragment->SetNumberOfParticles(excitons-holes);
  fragment->SetNumberOfParticles(excitons);
  fragment->SetNumberOfCharged(zCaptured);

#ifdef debug_G4BinaryCascade
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

  return fragment;
}

//----------------------------------------------------------------------------

G4LorentzVector G4BinaryCascade::GetFinal4Momentum()
//----------------------------------------------------------------------------
{
// the initial 3-momentum will differ from 0, if nucleus created by string model.
  G4LorentzVector final4Momentum = theInitial4Mom;
//G4cout << "GetFinal4Momentum:theInitial4Mom = " <<theInitial4Mom << G4endl;
//  G4double mass=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(the3DNucleus->GetCharge(), the3DNucleus->GetMassNumber());
//  G4LorentzVector final4Momentum(0,0,0,mass);
  G4KineticTrackVector::iterator i;
  for(i = theProjectileList.begin() ; i != theProjectileList.end(); ++i)
  {
    final4Momentum += (*i)->GetTrackingMomentum();
    //G4cerr << "Initial state: "<<(*i)->Get4Momentum()<<G4endl;
  }
//G4cout << "GetFinal4Momentum:theInitial4Mom+projectiles  = " <<final4Momentum << G4endl;

  for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
  {
    final4Momentum -= (*i)->Get4Momentum();
    // G4cerr <<"Final state: "<<(*i)->Get4Momentum()<<G4endl;
  }
//G4cout << "GetFinal4Momentum:theInitial4Mom+projectiles-finals  = " <<final4Momentum << G4endl;

//G4cerr << "what 0 "<< theProjectileList.size()<<" "<<theFinalState.size()<<G4endl;
//G4cerr << "what 1 - "<<theInitial4Mom<<" "<<final4Momentum<<G4endl;
  if((final4Momentum.vect()/final4Momentum.e()).mag()>1.0 && currentA > 0)
  {
    G4cerr << G4endl;
    G4cerr << "G4BinaryCascade::GetFinal4Momentum - Fatal"<<G4endl;
    G4KineticTrackVector::iterator i;
    G4cerr <<" Initial nucleus "<<theInitial4Mom<<G4endl;
    for(i = theProjectileList.begin() ; i != theProjectileList.end(); ++i)
    {
      G4cerr << " Initial state: "<<(*i)->Get4Momentum()<<", "<<(*i)->GetDefinition()->GetParticleName()<<G4endl;
    }
    for(i = theFinalState.begin(); i != theFinalState.end(); ++i)
    {
      G4cerr <<" Final state: "<<(*i)->Get4Momentum()<<(*i)->GetDefinition()->GetParticleName()<<G4endl;
    }
    G4cerr<< " Final4Momentum = "<<final4Momentum <<" "<<final4Momentum.m()<<G4endl;
    G4cerr <<" current A, Z = "<< currentA<<", "<<currentZ<<G4endl;
    G4cerr << G4endl;
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
// G4cout << "GetFinalNucleusMomentum GetFinal4Momentum= " <<NucleusMomentum <<" "<<NucleusMomentum.mag()<<G4endl;
// boost nucleus to a frame such that the momentum of nucleus == momentum of Captured
  G4ThreeVector boost= (NucleusMomentum.vect() -CapturedMomentum.vect())/NucleusMomentum.e();
//G4cout << "GetFinalNucleusMomentum boost= " <<boost<<G4endl;
#ifdef debug_G4BinaryCascade
  if(boost.mag()>1.0)
  {
    G4cerr << "G4BinaryCascade::GetFinalNucleusMomentum - Fatal"<<G4endl;
    G4cerr << "it 0"<<boost <<G4endl;
    G4cerr << "it 01"<<NucleusMomentum<<" "<<CapturedMomentum<<" "<<G4endl;
    G4cout <<" testing boost "<<boost<<" "<<boost.mag()<<G4endl;
  }
#endif
  G4LorentzRotation  nucleusBoost( -boost );
  precompoundLorentzboost.set( boost );
#ifdef debug_G4BinaryCascade
  G4cout << "it "<<NucleusMomentum<<" "<<CapturedMomentum<<" "<<G4endl;
 #endif
 NucleusMomentum *= nucleusBoost;
#ifdef debug_G4BinaryCascade
  G4cout << "GetFinalNucleusMomentum aft boost GetFinal4Momentum= " <<NucleusMomentum <<G4endl;
#endif
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
    while(!done)
    {
      if(secs)
      {
       std::for_each(secs->begin(), secs->end(), DeleteKineticTrack());
       delete secs;
      }
      secs = theScatterer->Scatter(*(*secondaries).front(), aTarget);
      for(size_t ss=0; secs && ss<secs->size(); ss++)
      {
        if((*secs)[ss]->GetDefinition()->IsShortLived()) done = true;
      }
    }
    size_t current(0);
    for(current=0; current<secs->size(); current++)
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
 *   G4double x = b*cos(phi);
 *   G4double y = b*sin(phi);
 *   G4double z = -sqrt(r*r-b*b);
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
    delete *i;
  ktv->clear();
}

//----------------------------------------------------------------------------
void G4BinaryCascade::ClearAndDestroy(G4ReactionProductVector * rpv)
//----------------------------------------------------------------------------
{
  std::vector<G4ReactionProduct *>::iterator i;
  for(i = rpv->begin(); i != rpv->end(); ++i)
    delete *i;
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
    G4cout << "  track n. " << count << ", id: " << kt << G4endl;
    G4ThreeVector pos = kt->GetPosition();
    G4LorentzVector mom = kt->Get4Momentum();
    G4ParticleDefinition * definition = kt->GetDefinition();
    G4cout << "    definition: " << definition->GetPDGEncoding() << " pos: "
	   << 1/fermi*pos << " R: " << 1/fermi*pos.mag() << " 4mom: "
	   << 1/MeV*mom << " P: " << 1/MeV*mom.vect().mag() 
	   << " M: " << 1/MeV*mom.mag() << G4endl;
    G4cout <<"trackstatus: "<<kt->GetState()<<G4endl;
  }
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
       G4LorentzVector pion(pion3, sqrt(sqr(140*MeV) +pion3.mag()));
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

void G4BinaryCascade::PrintWelcomeMessage()
{
  G4cout <<"Thank you for using G4BinaryCascade. "<<G4endl;
}








