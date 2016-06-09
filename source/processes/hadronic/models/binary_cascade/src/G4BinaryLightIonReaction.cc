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
#include "G4BinaryLightIonReaction.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include <algorithm>
#include "G4ReactionProductVector.hh"
#include <vector>
#include "G4ping.hh"
#include "G4Delete.hh"
#include "G4Neutron.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4HadTmpUtil.hh"
#include <cmath>
 
G4BinaryLightIonReaction::G4BinaryLightIonReaction()
    : G4HadronicInteraction("Binary Cascade"), theModel() , 
      theHandler(0) , theProjectileFragmentation(0)
{
    theHandler= new G4ExcitationHandler; 
    SetPrecompound(new G4PreCompoundModel(theHandler));
}
  
G4HadFinalState *G4BinaryLightIonReaction::
  ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus & targetNucleus )
{ 
  static G4int eventcounter=0;
  eventcounter++;
  if(getenv("BLICDEBUG") ) G4cerr << " ######### Binary Light Ion Reaction number starts ######### "<<eventcounter<<G4endl;
    G4ping debug("debug_G4BinaryLightIonReaction");
    G4int a1=aTrack.GetDefinition()->GetBaryonNumber();
    G4int z1=G4lrint(aTrack.GetDefinition()->GetPDGCharge());
    G4int a2=targetNucleus.GetA_asInt();
    G4int z2=targetNucleus.GetZ_asInt();
    debug.push_back(a1);
    debug.push_back(z1);
    debug.push_back(a2);
    debug.push_back(z2);
//    debug.push_back(m2);
    G4LorentzVector mom(aTrack.Get4Momentum());
    debug.push_back(mom);
    debug.dump();
    G4LorentzRotation toBreit(mom.boostVector());
        
    G4bool swapped = false;
    if(a2<a1)
    {
      swapped = true;
      G4int tmp(0);
      tmp = a2; a2=a1; a1=tmp;
      tmp = z2; z2=z1; z1=tmp;
      G4double m1=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z1,a1);
      G4LorentzVector it(m1, G4ThreeVector(0,0,0));
      mom = toBreit*it;
    }

    G4ReactionProductVector * result = NULL;
    G4ReactionProductVector * cascaders= new G4ReactionProductVector;
    G4double m_nucl(0);      // to check energy balance 


//    G4double m1=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z1,a1);
//    G4cout << "Entering the decision point "
//           << (mom.t()-mom.mag())/a1 << " "
//	   << a1<<" "<< z1<<" "
//	   << a2<<" "<< z2<<G4endl
//	   << " "<<mom.t()-mom.mag()<<" "
//	   << mom.t()- m1<<G4endl;
    if( (mom.t()-mom.mag())/a1 < 50*MeV )
    {
//      G4cout << "Using pre-compound only, E= "<<mom.t()-mom.mag()<<G4endl;
//      m_nucl = mom.mag();
      delete cascaders;
      G4Fragment aPreFrag;
      aPreFrag.SetA(a1+a2);
      aPreFrag.SetZ(z1+z2);
      aPreFrag.SetNumberOfParticles(a1);
      aPreFrag.SetNumberOfCharged(z1);
      aPreFrag.SetNumberOfHoles(0);
      G4ThreeVector plop(0.,0., mom.vect().mag());
      G4double m2=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z2,a2);
      m_nucl=m2;
      G4LorentzVector aL(mom.t()+m2, plop);
      aPreFrag.SetMomentum(aL);
      G4ParticleDefinition * preFragDef;
      preFragDef = G4ParticleTable::GetParticleTable()
                      ->FindIon(z1+z2,a1+a2,0,z1+z2);  
      aPreFrag.SetParticleDefinition(preFragDef);

//      G4cout << "Fragment INFO "<< a1+a2 <<" "<<z1+z2<<" "
//             << aL <<" "<<preFragDef->GetParticleName()<<G4endl;
      cascaders = theProjectileFragmentation->DeExcite(aPreFrag);
      G4double tSum = 0;
      for(size_t count = 0; count<cascaders->size(); count++)
      {
	cascaders->operator[](count)->SetNewlyAdded(true);
	tSum += cascaders->operator[](count)->GetKineticEnergy();
      }
//       G4cout << "Exiting pre-compound only, E= "<<tSum<<G4endl;
   }
    else
    {
 

      G4V3DNucleus * fancyNucleus = NULL;
      G4Fancy3DNucleus * projectile = NULL;
      G4double m1(0) ,m2(0);    
      G4LorentzVector it;

      G4FermiMomentum theFermi;
      G4int tryCount(0);
      while(!result)
      {
	projectile = new G4Fancy3DNucleus;
	projectile->Init(a1, z1);
	projectile->CenterNucleons();
	m1=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(
			  projectile->GetCharge(),projectile->GetMassNumber());
	it=toBreit * G4LorentzVector(m1,G4ThreeVector(0,0,0));
	fancyNucleus = new G4Fancy3DNucleus;  
	fancyNucleus->Init(a2, z2);
	m2=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(
      			  fancyNucleus->GetCharge(),fancyNucleus->GetMassNumber());
	m_nucl = ( swapped ) ? m1 : m2;
//	  G4cout << " mass table, nucleus, delta : " << m2 <<" "<< fancyNucleus->GetMass()
//               <<" "<<m2-fancyNucleus->GetMass() << G4endl;
	G4double impactMax = fancyNucleus->GetOuterRadius()+projectile->GetOuterRadius();
//        G4cout << "out radius - nucleus - projectile " << fancyNucleus->GetOuterRadius()/fermi << " - " << projectile->GetOuterRadius()/fermi << G4endl;
	G4double aX=(2.*G4UniformRand()-1.)*impactMax;
	G4double aY=(2.*G4UniformRand()-1.)*impactMax;
	G4ThreeVector pos(aX, aY, -2.*impactMax-5.*fermi);
	debug.push_back("Impact parameter");
	debug.push_back(aX);
	debug.push_back(aY);
	debug.push_back(-2.*impactMax);
	debug.dump();

	G4KineticTrackVector * initalState = new G4KineticTrackVector;
	projectile->StartLoop();
	G4Nucleon * aNuc;
	G4LorentzVector tmpV(0,0,0,0);
	G4LorentzVector nucleonMom(1./a1*mom);
	nucleonMom.setZ(nucleonMom.vect().mag());
	nucleonMom.setX(0);
	nucleonMom.setY(0);
	debug.push_back(" projectile nucleon momentum");
	debug.push_back(nucleonMom);
	debug.dump();
	theFermi.Init(a1,z1);
	while( (aNuc=projectile->GetNextNucleon()) )
	{
	  G4LorentzVector p4 = aNuc->GetMomentum();	
	  tmpV+=p4;
	  G4ThreeVector nucleonPosition(aNuc->GetPosition());
          G4double density=(projectile->GetNuclearDensity())->GetDensity(nucleonPosition);
	  nucleonPosition += pos;
	  G4KineticTrack * it = new G4KineticTrack(aNuc, nucleonPosition, nucleonMom );
          it->SetState(G4KineticTrack::outside);
	  G4double pfermi= theFermi.GetFermiMomentum(density);
	  G4double mass = aNuc->GetDefinition()->GetPDGMass();
	  G4double Efermi= std::sqrt( sqr(mass) + sqr(pfermi)) - mass;
          it->SetProjectilePotential(-Efermi);
	  initalState->push_back(it);
	}
	debug.push_back(" Sum of proj. nucleon momentum");
	debug.push_back(tmpV);
	debug.dump();

	result=theModel.Propagate(initalState, fancyNucleus);
	debug.push_back("################# Result size");
	if (result) {
	   debug.push_back(result->size());
	} else  {
	      debug.push_back(" -none-");
	}      
	debug.dump();
//	std::for_each(initalState->begin(), initalState->end(), Delete<G4KineticTrack>());
//	delete initalState;

	if(! result || result->size()==0) 
	{
          if (result) {delete result; result=0;}
          delete fancyNucleus;
          delete projectile;
	  if (++tryCount > 200)
	  {
	      // abort!!
	      
	      G4cerr << "G4BinaryLightIonReaction no final state for: " << G4endl;
	      G4cerr << " Primary " << aTrack.GetDefinition()
	  	       << ", (A,Z)=(" << aTrack.GetDefinition()->GetBaryonNumber()
		       << "," << aTrack.GetDefinition()->GetPDGCharge() << ") "
	               << ", kinetic energy " << aTrack.GetKineticEnergy() 
		       << G4endl;
	      G4cerr << " Target nucleus (A,Z)=(" <<  targetNucleus.GetA_asInt()
	               << "," << targetNucleus.GetZ_asInt() << G4endl;
	      G4cerr << " if frequent, please submit above information as bug report"
  		      << G4endl << G4endl;
		
	      theResult.Clear();
	      theResult.SetStatusChange(isAlive);
	      theResult.SetEnergyChange(aTrack.GetKineticEnergy());
	      theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    	      return &theResult;

	  }
	} 
	else
	{
          break;
	}     
      }
	debug.push_back(" Attempts to create final state");
	debug.push_back(tryCount);
	debug.dump();
      debug.push_back("################# Through the loop ? "); debug.dump();
      //inverse transformation in case we swapped.
      G4int resA(0), resZ(0); 
      G4Nucleon * aNuc;
//       fancyNucleus->StartLoop();
//       while( (aNuc=fancyNucleus->GetNextNucleon()) )
//       {
//         G4cout << " tgt Nucleon : " << aNuc->GetDefinition()->GetParticleName() <<" "<< aNuc->AreYouHit() <<" "<<aNuc->GetMomentum()<<G4endl;
//       }
      G4ReactionProductVector * spectators= new G4ReactionProductVector;
     debug.push_back("getting at the hits"); debug.dump();
      // the projectile excitation energy estimate...
      G4double theStatisticalExEnergy = 0;
      projectile->StartLoop();  
      while( (aNuc=projectile->GetNextNucleon()) )
      {
//        G4cout << " Nucleon : " << aNuc->GetDefinition()->GetParticleName() <<" "<< aNuc->AreYouHit() <<" "<<aNuc->GetMomentum()<<G4endl;
	debug.push_back("getting the hits"); debug.dump();
	if(!aNuc->AreYouHit())
	{
          resA++;
	  resZ+=G4lrint(aNuc->GetDefinition()->GetPDGCharge());
	}
	else
	{
          debug.push_back(" ##### a hit ##### "); debug.dump();
	  G4ThreeVector aPosition(aNuc->GetPosition());
          G4double localDensity = projectile->GetNuclearDensity()->GetDensity(aPosition);
	  G4double localPfermi = theFermi.GetFermiMomentum(localDensity);
	  G4double nucMass = aNuc->GetDefinition()->GetPDGMass();
	  G4double localFermiEnergy = std::sqrt(nucMass*nucMass + localPfermi*localPfermi) - nucMass;
	  G4double deltaE = localFermiEnergy - (aNuc->GetMomentum().t()-aNuc->GetMomentum().mag());
	  theStatisticalExEnergy += deltaE;
	}
	debug.push_back("collected a hit"); 
	debug.push_back(aNuc->GetMomentum());
	debug.dump();
      }
      delete fancyNucleus;
      delete projectile;
      G4ping debug("debug_G4BinaryLightIonReaction_1");
      debug.push_back("have the hits. A,Z, excitE"); 
      debug.push_back(resA);
      debug.push_back(resZ);
      debug.push_back(theStatisticalExEnergy);
      debug.dump();
      // Calculate excitation energy
      G4LorentzVector iState = mom;
      iState.setT(iState.getT()+m2);

      G4LorentzVector fState(0,0,0,0);
      G4LorentzVector pspectators(0,0,0,0);
      unsigned int i(0);
//      G4int spectA(0),spectZ(0);
      for(i=0; i<result->size(); i++)
      {
	if( (*result)[i]->GetNewlyAdded() ) 
	{
          fState += G4LorentzVector( (*result)[i]->GetMomentum(), (*result)[i]->GetTotalEnergy() );
	  cascaders->push_back((*result)[i]);
//          G4cout <<" secondary ... ";
          debug.push_back("secondary ");
	  debug.push_back((*result)[i]->GetDefinition()->GetParticleName());
	  debug.push_back(G4LorentzVector((*result)[i]->GetMomentum(),(*result)[i]->GetTotalEnergy()));
	  debug.dump();
	}
	else {
//          G4cout <<" spectator ... ";
          pspectators += G4LorentzVector( (*result)[i]->GetMomentum(), (*result)[i]->GetTotalEnergy() );
	  spectators->push_back((*result)[i]);
          debug.push_back("spectator ");
	  debug.push_back((*result)[i]->GetDefinition()->GetParticleName());
	  debug.push_back(G4LorentzVector((*result)[i]->GetMomentum(),(*result)[i]->GetTotalEnergy()));
	  debug.dump();
//	  spectA++; 
//	  spectZ+= G4lrint((*result)[i]->GetDefinition()->GetPDGCharge());
	}

//       G4cout << (*result)[i]<< " "
//   	    << (*result)[i]->GetDefinition()->GetParticleName() << " " 
//   	    << (*result)[i]->GetMomentum()<< " " 
//   	    << (*result)[i]->GetTotalEnergy() << G4endl;
      }
//      if ( spectA-resA !=0 || spectZ-resZ !=0)
//      {
//          G4cout << "spect Nucl != spectators: nucl a,z; spect a,z" <<
//	      resA <<" "<< resZ <<" ; " << spectA <<" "<< spectZ << G4endl;
//      }
      delete result;
      debug.push_back(" iState - (fState+pspectators) ");
      debug.push_back(iState-fState-pspectators);
      debug.dump();
      G4LorentzVector momentum(iState-fState);
      G4int loopcount(0);
      while (std::abs(momentum.e()-pspectators.e()) > 10*MeV)
      {
	 debug.push_back("the momentum balance");
	 debug.push_back(iState);
	 debug.push_back(fState);
	 debug.push_back(momentum-pspectators);
	 debug.push_back(momentum);
	 debug.dump();
	 G4LorentzVector pCorrect(iState-pspectators);
	 G4bool EnergyIsCorrect=EnergyAndMomentumCorrector(cascaders, pCorrect);
         if ( ! EnergyIsCorrect && getenv("debug_G4BinaryLightIonReactionResults"))
         {
            G4cout << "Warning - G4BinaryLightIonReaction E/P correction for cascaders failed" << G4endl;
         }
	 fState=G4LorentzVector();
	 for(i=0; i<cascaders->size(); i++)
         {
	     fState += G4LorentzVector( (*cascaders)[i]->GetMomentum(), (*cascaders)[i]->GetTotalEnergy() );
	 }
	 momentum=iState-fState;
	 debug.push_back("the momentum balance after correction");
	 debug.push_back(iState);
	 debug.push_back(fState);
	 debug.push_back(momentum-pspectators);
	 debug.push_back(momentum);
	 debug.dump();
	 if (++loopcount > 10 ) 
	 {   
	     if ( momentum.vect().mag() > momentum.e() )
	     {
	        G4cerr << "G4BinaryLightIonReaction.cc: Cannot correct 4-momentum of cascade particles" << G4endl;
	        throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde::ApplyCollision()");
	     } else {
	        break;
	     }
	     
 	 }
      }




      // call precompound model
      G4ReactionProductVector * proFrag = NULL;
      G4LorentzVector pFragment;
//      G4cout << " == pre boost 1 "<< momentum.e()<< " "<< momentum.mag()<<G4endl;
      G4LorentzRotation boost_fragments;
//      G4cout << " == post boost 1 "<< momentum.e()<< " "<< momentum.mag()<<G4endl;
  //    G4LorentzRotation boost_spectator_mom(-momentum.boostVector());
  //     G4cout << "- momentum " << boost_spectator_mom * momentum << G4endl; 
      G4LorentzVector pFragments(0);
      if(resZ>0 && resA>1) 
      {
	//  Make the fragment
	G4Fragment aProRes;
	aProRes.SetA(resA);
	aProRes.SetZ(resZ);
	aProRes.SetNumberOfParticles(0);
	aProRes.SetNumberOfCharged(0);
	aProRes.SetNumberOfHoles(a1-resA);
	G4double mFragment=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(resZ,resA);
	G4LorentzVector pFragment(0,0,0,mFragment+std::max(0.,theStatisticalExEnergy) );
	aProRes.SetMomentum(pFragment);
	G4ParticleDefinition * resDef;
	resDef = G4ParticleTable::GetParticleTable()->FindIon(resZ,resA,0,resZ);  
	aProRes.SetParticleDefinition(resDef);
	proFrag = theHandler->BreakItUp(aProRes);
      if ( momentum.vect().mag() > momentum.e() )
      {
           G4cerr << "mom check: " <<  momentum 
	          << " 3.mag "<< momentum.vect().mag() << G4endl
		  << " .. iState/fState/spectators " << iState <<" " 
		  << fState << " " << pspectators << G4endl
		  << " .. A,Z " << resA <<" "<< resZ << G4endl;          	                
	   G4cerr << "G4BinaryLightIonReaction no final state for: " << G4endl;
	   G4cerr << " Primary " << aTrack.GetDefinition()
	  	    << ", (A,Z)=(" << aTrack.GetDefinition()->GetBaryonNumber()
		    << "," << aTrack.GetDefinition()->GetPDGCharge() << ") "
	            << ", kinetic energy " << aTrack.GetKineticEnergy() 
		    << G4endl;
	   G4cerr << " Target nucleus (A,Z)=(" <<  targetNucleus.GetA_asInt()
	            << "," << targetNucleus.GetZ_asInt() << G4endl;
	   G4cerr << " if frequent, please submit above information as bug report"
  	     	   << G4endl << G4endl;

	   theResult.Clear();
	   theResult.SetStatusChange(isAlive);
	   theResult.SetEnergyChange(aTrack.GetKineticEnergy());
	   theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
    	   return &theResult;
      }
      
        G4LorentzRotation boost_fragments_here(momentum.boostVector());
	boost_fragments = boost_fragments_here;

//   	 G4cout << " Fragment a,z, Mass Fragment, mass spect-mom, exitationE " 
//   		<< resA <<" "<< resZ <<" "<< mFragment <<" "
//   		<< momentum.mag() <<" "<< momentum.mag() - mFragment 
//   	   << " "<<theStatisticalExEnergy 
// 	   << " "<< boost_fragments*pFragment<< G4endl;
        G4ReactionProductVector::iterator ispectator;
	for (ispectator=spectators->begin();ispectator!=spectators->end();ispectator++)
	{
	    delete *ispectator;
	}
      }
      else if(resA!=0)
      {
           G4ReactionProductVector::iterator ispectator;
	   for (ispectator=spectators->begin();ispectator!=spectators->end();ispectator++)
	   {
	       (*ispectator)->SetNewlyAdded(true);
	       cascaders->push_back(*ispectator);
	       pFragments+=G4LorentzVector((*ispectator)->GetMomentum(),(*ispectator)->GetTotalEnergy());
  // 	     G4cout << "from spectator "  
  //  	      << (*ispectator)->GetDefinition()->GetParticleName() << " " 
  //  	      << (*ispectator)->GetMomentum()<< " " 
  //  	      << (*ispectator)->GetTotalEnergy() << G4endl;
	   }
      }
      if (spectators) delete spectators;

      // collect the evaporation part
      debug.push_back("the nucleon count balance");
      debug.push_back(resA);
      debug.push_back(resZ);
      if(proFrag) debug.push_back(proFrag->size());
      debug.dump();
      G4ReactionProductVector::iterator ii;
      if(proFrag) for(ii=proFrag->begin(); ii!=proFrag->end(); ii++)
      {
	(*ii)->SetNewlyAdded(true);
	G4LorentzVector tmp((*ii)->GetMomentum(),(*ii)->GetTotalEnergy());
	tmp *= boost_fragments;
	(*ii)->SetMomentum(tmp.vect());
	(*ii)->SetTotalEnergy(tmp.e());
  //      result->push_back(*ii);
	pFragments += tmp;
      }

  //    G4cout << "Fragmented p, momentum, delta " << pFragments <<" "<<momentum
  //            <<" "<< pFragments-momentum << G4endl;
      debug.push_back("################# done with evaporation"); debug.dump();

  //  correct p/E of Cascade secondaries
      G4LorentzVector pCas=iState - pFragments;
  //    G4cout <<" Going to correct from " << fState << " to " << pCas << G4endl;
      G4bool EnergyIsCorrect=EnergyAndMomentumCorrector(cascaders, pCas);
      if ( ! EnergyIsCorrect )
      {
	 if(getenv("debug_G4BinaryLightIonReactionResults"))      
           G4cout << "G4BinaryLightIonReaction E/P correction for nucleus failed, will try to correct overall" << G4endl;
      }

  //  Add deexcitation secondaries 
      if(proFrag) for(ii=proFrag->begin(); ii!=proFrag->end(); ii++)
      {
	cascaders->push_back(*ii);
      }
      if (proFrag) delete proFrag;

      if ( ! EnergyIsCorrect )
      {
         if (! EnergyAndMomentumCorrector(cascaders,iState))
	 { 
	   if(getenv("debug_G4BinaryLightIonReactionResults"))      
             G4cout << "G4BinaryLightIonReaction E/P corrections failed" << G4endl;
	 }
      }
      
    }
    
    // Rotate to lab
    G4LorentzRotation toZ;
    toZ.rotateZ(-1*mom.phi());
    toZ.rotateY(-1*mom.theta());
    G4LorentzRotation toLab(toZ.inverse());
  
    // Fill the particle change, while rotating. Boost from projectile breit-frame in case we swapped.  
    // theResult.Clear();
    theResult.Clear();
    theResult.SetStatusChange(stopAndKill);
    G4double Etot(0);
    size_t i=0;
    for(i=0; i<cascaders->size(); i++)
    {
      if((*cascaders)[i]->GetNewlyAdded())
      {
	G4DynamicParticle * aNew = 
	new G4DynamicParticle((*cascaders)[i]->GetDefinition(),
                              (*cascaders)[i]->GetTotalEnergy(),
			      (*cascaders)[i]->GetMomentum() );
	G4LorentzVector tmp = aNew->Get4Momentum();
	if(swapped)
	{
          tmp*=toBreit.inverse();
	  tmp.setVect(-tmp.vect());
	}    
	tmp *= toLab;
	aNew->Set4Momentum(tmp);
	theResult.AddSecondary(aNew);
	Etot += tmp.e();
//        G4cout << "LIBIC: Secondary " << aNew->GetDefinition()->GetParticleName()
//               <<" "<<  aNew->GetMomentum()
// 	      <<" "<<  aNew->GetTotalEnergy()
// 	      << G4endl;
      }
      delete (*cascaders)[i];
    }
    if(cascaders) delete cascaders;
    
    G4ping debug1("debug_G4BinaryLightIonReactionResults");
    debug1.push_back("Result analysis, secondaries");
    debug1.push_back(theResult.GetNumberOfSecondaries());
    debug1.dump();
    debug1.push_back(" Energy conservation initial/final/delta(init-final) ");
    debug1.push_back(aTrack.GetTotalEnergy() + m_nucl);
    debug1.push_back(aTrack.GetTotalEnergy());
    debug1.push_back(m_nucl);
    debug1.push_back(Etot);
    debug1.push_back(aTrack.GetTotalEnergy() + m_nucl - Etot);
    debug1.dump();

    if(getenv("BLICDEBUG") ) G4cerr << " ######### Binary Light Ion Reaction number ends ######### "<<eventcounter<<G4endl;

    return &theResult;
  }

//****************************************************************************  
G4bool G4BinaryLightIonReaction::EnergyAndMomentumCorrector(
	G4ReactionProductVector* Output, G4LorentzVector& TotalCollisionMom)   
//****************************************************************************  
  {
    const int    nAttemptScale = 2500;
    const double ErrLimit = 1.E-6;
    if (Output->empty())
       return TRUE;
    G4LorentzVector SumMom(0);
    G4double        SumMass = 0;     
    G4double        TotalCollisionMass = TotalCollisionMom.m(); 
    size_t i = 0;
    // Calculate sum hadron 4-momenta and summing hadron mass
    for(i = 0; i < Output->size(); i++)
        {
        SumMom  += G4LorentzVector((*Output)[i]->GetMomentum(),(*Output)[i]->GetTotalEnergy());
        SumMass += (*Output)[i]->GetDefinition()->GetPDGMass();
        }
//   G4cout << " E/P corrector, SumMass, SumMom.m2, TotalMass " 
//         << SumMass <<" "<< SumMom.m2() <<" "<<TotalCollisionMass<< G4endl; 
    if (SumMass > TotalCollisionMass) return FALSE;
    SumMass = SumMom.m2();
    if (SumMass < 0) return FALSE;
    SumMass = std::sqrt(SumMass);

     // Compute c.m.s. hadron velocity and boost KTV to hadron c.m.s.
    G4ThreeVector Beta = -SumMom.boostVector();
//      G4cout << " == pre boost 2 "<< SumMom.e()<< " "<< SumMom.mag()<<" "<< Beta <<G4endl;
//--old    Output->Boost(Beta);
      for(i = 0; i < Output->size(); i++)
      {
        G4LorentzVector mom = G4LorentzVector((*Output)[i]->GetMomentum(),(*Output)[i]->GetTotalEnergy());
        mom *= Beta;
        (*Output)[i]->SetMomentum(mom.vect());
        (*Output)[i]->SetTotalEnergy(mom.e());
      }   

    // Scale total c.m.s. hadron energy (hadron system mass).
    // It should be equal interaction mass
    G4double Scale = 0,OldScale=0;
    G4double factor = 1.;
    G4int cAttempt = 0;
    G4double Sum = 0;
    G4bool success = false;
    for(cAttempt = 0; cAttempt < nAttemptScale; cAttempt++)
    {
      Sum = 0;
      for(i = 0; i < Output->size(); i++)
      {
        G4LorentzVector HadronMom = G4LorentzVector((*Output)[i]->GetMomentum(),(*Output)[i]->GetTotalEnergy());
        HadronMom.setVect(HadronMom.vect()+ factor*Scale*HadronMom.vect());
        G4double E = std::sqrt(HadronMom.vect().mag2() + sqr((*Output)[i]->GetDefinition()->GetPDGMass()));
        HadronMom.setE(E);
        (*Output)[i]->SetMomentum(HadronMom.vect());
        (*Output)[i]->SetTotalEnergy(HadronMom.e());
        Sum += E;
      }
      OldScale=Scale;   
      Scale = TotalCollisionMass/Sum - 1;
//  G4cout << "E/P corr - " << cAttempt << " " << Scale << G4endl;
      if (std::abs(Scale) <= ErrLimit  
          || OldScale == Scale)			// protect 'frozen' situation and divide by 0 in calculating new factor below
      {
        if (getenv("debug_G4BinaryLightIonReactionResults")) G4cout << "E/p corrector: " << cAttempt << G4endl;
        success = true;
	break;
      }
      if ( cAttempt > 10 ) 
      {
//         G4cout << " speed it up? " << std::abs(OldScale/(OldScale-Scale)) << G4endl;
	 factor=std::max(1.,std::log(std::abs(OldScale/(OldScale-Scale))));
//	 G4cout << " ? factor ? " << factor << G4endl;
      }   
    }
    
    if( (!success)  && getenv("debug_G4BinaryLightIonReactionResults"))      
    {
      G4cout << "G4G4BinaryLightIonReaction::EnergyAndMomentumCorrector - Warning"<<G4endl;
      G4cout << "   Scale not unity at end of iteration loop: "<<TotalCollisionMass<<" "<<Sum<<" "<<Scale<<G4endl;
      G4cout << "   Increase number of attempts or increase ERRLIMIT"<<G4endl;
    }

    // Compute c.m.s. interaction velocity and KTV back boost   
    Beta = TotalCollisionMom.boostVector(); 
//--old    Output->Boost(Beta);
      for(i = 0; i < Output->size(); i++)
      {
        G4LorentzVector mom = G4LorentzVector((*Output)[i]->GetMomentum(),(*Output)[i]->GetTotalEnergy());
        mom *= Beta;
        (*Output)[i]->SetMomentum(mom.vect());
        (*Output)[i]->SetTotalEnergy(mom.e());
      }   
    return TRUE;
  }
