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
#include <algorithm>
#include <vector>
#include <cmath>
#include <numeric>

#include "G4BinaryLightIonReaction.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include "G4ReactionProductVector.hh"
#include "G4ping.hh"
#include "G4Delete.hh"
#include "G4Neutron.hh"
#include "G4VNuclearDensity.hh"
#include "G4FermiMomentum.hh"
#include "G4HadTmpUtil.hh"
#include "G4PreCompoundModel.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4Log.hh"

#ifdef G4MULTITHREADED
   G4Mutex G4BinaryLightIonReaction::BLIRMutex = G4MUTEX_INITIALIZER;
#endif
   G4int G4BinaryLightIonReaction::theBLIR_ID = -1;


//#define debug_G4BinaryLightIonReaction
//#define debug_BLIR_finalstate
//#define debug_BLIR_result

G4BinaryLightIonReaction::G4BinaryLightIonReaction(G4VPreCompoundModel* ptr)
: G4HadronicInteraction("Binary Light Ion Cascade"),
  theProjectileFragmentation(ptr),
  pA(0),pZ(0), tA(0),tZ(0),spectatorA(0),spectatorZ(0),
  projectile3dNucleus(0),target3dNucleus(0)
{
	if(!ptr) {
	  G4HadronicInteraction* p =
	    G4HadronicInteractionRegistry::Instance()->FindModel("PRECO");
	  G4VPreCompoundModel* pre = static_cast<G4VPreCompoundModel*>(p);
	  if(!pre) { pre = new G4PreCompoundModel(); }
	  theProjectileFragmentation = pre;
	}
	theModel = new G4BinaryCascade(theProjectileFragmentation);
	theHandler = theProjectileFragmentation->GetExcitationHandler();
	if ( theBLIR_ID == -1 ) {
#ifdef G4MULTITHREADED
       G4MUTEXLOCK(&G4BinaryLightIonReaction::BLIRMutex);
       if ( theBLIR_ID == -1 ) {
#endif
    	   theBLIR_ID = G4PhysicsModelCatalog::Register("Binary Light Ion Reaction");
#ifdef G4MULTITHREADED
       }
       G4MUTEXUNLOCK(&G4BinaryLightIonReaction::BLIRMutex);
#endif
    }


	debug_G4BinaryLightIonReactionResults=getenv("debug_G4BinaryLightIonReactionResults")!=0;
}

G4BinaryLightIonReaction::~G4BinaryLightIonReaction()
{}

void G4BinaryLightIonReaction::ModelDescription(std::ostream& outFile) const
{
	outFile << "G4Binary Light Ion Cascade is an intra-nuclear cascade model\n"
			<< "using G4BinaryCasacde to model the interaction of a light\n"
			<< "nucleus with a nucleus.\n"
			<< "The lighter of the two nuclei is treated like a set of projectiles\n"
			<< "which are transported simultanously through the heavier nucleus.\n";
}

//--------------------------------------------------------------------------------
struct ReactionProduct4Mom
{
   G4LorentzVector operator()(G4LorentzVector a,G4ReactionProduct* b) {return a + G4LorentzVector(b->GetMomentum(), b->GetTotalEnergy() );}
};

G4HadFinalState *G4BinaryLightIonReaction::
ApplyYourself(const G4HadProjectile &aTrack, G4Nucleus & targetNucleus )
{
	if(getenv("BLICDEBUG") ) G4cerr << " ######### Binary Light Ion Reaction starts ######### " << G4endl;
	G4ping debug("debug_G4BinaryLightIonReaction");
	pA=aTrack.GetDefinition()->GetBaryonNumber();
	pZ=G4lrint(aTrack.GetDefinition()->GetPDGCharge()/eplus);
	tA=targetNucleus.GetA_asInt();
	tZ=targetNucleus.GetZ_asInt();
	G4double timePrimary = aTrack.GetGlobalTime();
	G4LorentzVector mom(aTrack.Get4Momentum());
   //G4cout << "proj mom : " << mom << G4endl;
	G4LorentzRotation toBreit(mom.boostVector());

	G4bool swapped=SetLighterAsProjectile(mom, toBreit);
   //G4cout << "after swap, swapped? / mom " << swapped << " / " << mom <<G4endl;
	G4ReactionProductVector * result = 0;
	G4ReactionProductVector * cascaders=0; //new G4ReactionProductVector;
//	G4double m_nucl(0);      // to check energy balance


	//    G4double m1=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(pZ,pA);
	//    G4cout << "Entering the decision point "
	//           << (mom.t()-mom.mag())/pA << " "
	//	   << pA<<" "<< pZ<<" "
	//	   << tA<<" "<< tZ<<G4endl
	//	   << " "<<mom.t()-mom.mag()<<" "
	//	   << mom.t()- m1<<G4endl;
	if( (mom.t()-mom.mag())/pA < 50*MeV )
	{
		//      G4cout << "Using pre-compound only, E= "<<mom.t()-mom.mag()<<G4endl;
		//      m_nucl = mom.mag();
      cascaders=FuseNucleiAndPrompound(mom);
      if( !cascaders )
      {

         // abort!! happens for too low energy for nuclei to fuse

         theResult.Clear();
         theResult.SetStatusChange(isAlive);
         theResult.SetEnergyChange(aTrack.GetKineticEnergy());
         theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
         return &theResult;
      }
	}
	else
	{
	   result=Interact(mom,toBreit);

      if(! result )
      {
            // abort!!

            G4cerr << "G4BinaryLightIonReaction no final state for: " << G4endl;
            G4cerr << " Primary " << aTrack.GetDefinition()
               << ", (A,Z)=(" << aTrack.GetDefinition()->GetBaryonNumber()
               << "," << aTrack.GetDefinition()->GetPDGCharge()/eplus << ") "
               << ", kinetic energy " << aTrack.GetKineticEnergy()
               << G4endl;
            G4cerr << " Target nucleus (A,Z)=("
                   <<  (swapped?pA:tA)  << ","
                   << (swapped?pZ:tZ) << ")" << G4endl;
            G4cerr << " if frequent, please submit above information as bug report"
                  << G4endl << G4endl;

            theResult.Clear();
            theResult.SetStatusChange(isAlive);
            theResult.SetEnergyChange(aTrack.GetKineticEnergy());
            theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
            return &theResult;
      }

		// Calculate excitation energy,
      G4double theStatisticalExEnergy = GetProjectileExcitation();


		pInitialState = mom;
		   //G4cout << "BLIC: pInitialState from aTrack : " << pInitialState;
		pInitialState.setT(pInitialState.getT() +
		     G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(tZ,tA));
           //G4cout << "BLIC: target nucleus added : " << pInitialState << G4endl;

		delete target3dNucleus;target3dNucleus=0;
		delete projectile3dNucleus;projectile3dNucleus=0;

      G4ReactionProductVector * spectators= new G4ReactionProductVector;

      cascaders = new G4ReactionProductVector;

      G4LorentzVector pspectators=SortResult(result,spectators,cascaders);
        // this also sets spectatorA and spectatorZ

      //      pFinalState=std::accumulate(cascaders->begin(),cascaders->end(),pFinalState,ReactionProduct4Mom);

      std::vector<G4ReactionProduct *>::iterator iter;

        // G4cout << "pInitialState, pFinalState / pspectators"<< pInitialState << " / " << pFinalState << " / " << pspectators << G4endl;
		//      if ( spectA-spectatorA !=0 || spectZ-spectatorZ !=0)
		//      {
		//          G4cout << "spect Nucl != spectators: nucl a,z; spect a,z" <<
		//	      spectatorA <<" "<< spectatorZ <<" ; " << spectA <<" "<< spectZ << G4endl;
		//      }
		delete result;
		result=0;
		G4LorentzVector momentum(pInitialState-pFinalState);
		G4int loopcount(0);
		   //G4cout << "BLIC: momentum, pspectators : " << momentum << " / " << pspectators << G4endl;
		while (std::abs(momentum.e()-pspectators.e()) > 10*MeV)                /* Loop checking, 31.08.2015, G.Folger */
		                                                                       // see if on loopcount 
		{
			G4LorentzVector pCorrect(pInitialState-pspectators);
		         //G4cout << "BLIC:: BIC nonconservation? (pInitialState-pFinalState) / spectators :" << momentum << " / " << pspectators << "pCorrect "<< pCorrect<< G4endl;
			// Correct outgoing casacde particles.... to have momentum of (initial state - spectators)
			G4bool EnergyIsCorrect=EnergyAndMomentumCorrector(cascaders, pCorrect);
			if ( ! EnergyIsCorrect && debug_G4BinaryLightIonReactionResults)
			{
				G4cout << "Warning - G4BinaryLightIonReaction E/P correction for cascaders failed" << G4endl;
			}
			pFinalState=G4LorentzVector(0,0,0,0);
		  	for(iter=cascaders->begin(); iter!=cascaders->end(); iter++)
			{
				pFinalState += G4LorentzVector( (*iter)->GetMomentum(), (*iter)->GetTotalEnergy() );
			}
			momentum=pInitialState-pFinalState;
			if (++loopcount > 10 )
			{
				if ( momentum.vect().mag() - momentum.e()> 10*keV  )
				{
					G4cerr << "G4BinaryLightIonReaction.cc: Cannot correct 4-momentum of cascade particles" << G4endl;
					throw G4HadronicException(__FILE__, __LINE__, "G4BinaryCasacde::ApplyCollision()");
				} else {
					break;
				}

			}
		}

   	if (spectatorA > 0 )
		{
		   // check spectator momentum
		   if ( momentum.vect().mag() - momentum.e()< 10*keV )
		   {
   	   	   	   // DeExciteSpectatorNucleus() also handles also case of A=1, Z=0,1
               DeExciteSpectatorNucleus(spectators, cascaders, theStatisticalExEnergy, momentum);

		   } else {   // momentum non-conservation --> fail
			  for (iter=spectators->begin();iter!=spectators->end();iter++)
		      {
		    	 delete *iter;
		      }
		      delete spectators;
		  	  for(iter=cascaders->begin(); iter!=cascaders->end(); iter++)
		  	  {
		  		 delete *iter;
		  	  }
		  	  delete cascaders;

		      G4cout << "G4BinaryLightIonReaction.cc: mom check: " <<  momentum
		            << " 3.mag "<< momentum.vect().mag() << G4endl
		            << " .. pInitialState/pFinalState/spectators " << pInitialState <<" "
		            << pFinalState << " " << pspectators << G4endl
		            << " .. A,Z " << spectatorA <<" "<< spectatorZ << G4endl;
		      G4cout << "G4BinaryLightIonReaction invalid final state for: " << G4endl;
		      G4cout << " Primary " << aTrack.GetDefinition()
  	    		      << ", (A,Z)=(" << aTrack.GetDefinition()->GetBaryonNumber()
  	    		      << "," << aTrack.GetDefinition()->GetPDGCharge()/eplus << ") "
  	    		      << ", kinetic energy " << aTrack.GetKineticEnergy()
  	    		      << G4endl;
		      G4cout << " Target nucleus (A,Z)=(" <<  targetNucleus.GetA_asInt()
            		      << "," << targetNucleus.GetZ_asInt() << ")" << G4endl;
		      G4cout << " if frequent, please submit above information as bug report"
		            << G4endl << G4endl;
#ifdef debug_G4BinaryLightIonReaction
		      G4ExceptionDescription ed;
		      ed << "G4BinaryLightIonreaction: Terminate for above error"  << G4endl;
		      G4Exception("G4BinaryLightIonreaction::ApplyYourSelf()", "BLIC001", FatalException,
		  		ed);

#endif
		      theResult.Clear();
		      theResult.SetStatusChange(isAlive);
		      theResult.SetEnergyChange(aTrack.GetKineticEnergy());
		      theResult.SetMomentumChange(aTrack.Get4Momentum().vect().unit());
		      return &theResult;
		   }
		} else {    // no spectators
		   delete spectators;
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
	G4LorentzVector ptot(0);
    G4ReactionProductVector::iterator iter;
    #ifdef debug_BLIR_result
       G4LorentzVector p_raw;
    #endif
    //G4int i=0;

	for(iter=cascaders->begin(); iter!=cascaders->end(); iter++)
	{
		if((*iter)->GetNewlyAdded())
		{
			G4DynamicParticle * aNewDP =
					new G4DynamicParticle((*iter)->GetDefinition(),
							(*iter)->GetTotalEnergy(),
							(*iter)->GetMomentum() );
			G4LorentzVector tmp = aNewDP->Get4Momentum();
             #ifdef debug_BLIR_result
			     p_raw+= tmp;
             #endif
			if(swapped)
			{
				tmp*=toBreit.inverse();
				tmp.setVect(-tmp.vect());
			}
			tmp *= toLab;
			aNewDP->Set4Momentum(tmp);
			G4HadSecondary aNew = G4HadSecondary(aNewDP);
            G4double time = 0;                     //(*iter)->GetCreationTime();
            //if(time < 0.0) { time = 0.0; }
            aNew.SetTime(timePrimary + time);
            aNew.SetCreatorModelType((*iter)->GetCreatorModel());

			theResult.AddSecondary(aNew);
			ptot += tmp;
			        //G4cout << "BLIC: Secondary " << aNew->GetDefinition()->GetParticleName()
			        //       <<" "<<  aNew->GetMomentum()<<" "<<  aNew->GetTotalEnergy() << G4endl;
		}
		delete *iter;
	}
	delete cascaders;

#ifdef debug_BLIR_result
	//G4cout << "Result analysis, secondaries " << theResult.GetNumberOfSecondaries() << G4endl;
	//G4cout << "p_tot_raw " << p_raw << " sum p final " << ptot << G4endl;
	G4double m_nucl=	G4ParticleTable::GetParticleTable()->GetIonTable()->
			    GetIonMass(targetNucleus.GetZ_asInt(),targetNucleus.GetA_asInt());
	// delete? tZ=targetNucleus.GetZ_asInt();

	//G4cout << "BLIC Energy conservation initial/primary/nucleus/final/delta(init-final) "
	 //     << aTrack.GetTotalEnergy()   + m_nucl <<" "<< aTrack.GetTotalEnergy() <<" "<< m_nucl <<" "<<ptot.e()
	 //     <<" "<< aTrack.GetTotalEnergy() + m_nucl - ptot.e() << G4endl;
	G4cout << "BLIC momentum conservation " << aTrack.Get4Momentum()+ G4LorentzVector(m_nucl)
			<< " ptot " << ptot << " delta " << aTrack.Get4Momentum()+ G4LorentzVector(m_nucl) - ptot
			<< "        3mom.mag() " << (aTrack.Get4Momentum()+ G4LorentzVector(m_nucl) - ptot).vect().mag() << G4endl;
#endif

	if(getenv("BLICDEBUG") ) G4cerr << " ######### Binary Light Ion Reaction number ends ######### " << G4endl;

	return &theResult;
}

//--------------------------------------------------------------------------------

//****************************************************************************
G4bool G4BinaryLightIonReaction::EnergyAndMomentumCorrector(
		G4ReactionProductVector* Output, G4LorentzVector& TotalCollisionMom)
//****************************************************************************
{
	const int    nAttemptScale = 2500;
	const double ErrLimit = 1.E-6;
	if (Output->empty())
		return TRUE;
	G4LorentzVector SumMom(0,0,0,0);
	G4double        SumMass = 0;
	G4double        TotalCollisionMass = TotalCollisionMom.m();
	size_t i = 0;
	// Calculate sum hadron 4-momenta and summing hadron mass
	for(i = 0; i < Output->size(); i++)
	{
		SumMom  += G4LorentzVector((*Output)[i]->GetMomentum(),(*Output)[i]->GetTotalEnergy());
		SumMass += (*Output)[i]->GetDefinition()->GetPDGMass();
	}
	  // G4cout << " E/P corrector, SumMass, SumMom.m2, TotalMass "
	  //       << SumMass <<" "<< SumMom.m2() <<" "<<TotalCollisionMass<< G4endl;
	if (SumMass > TotalCollisionMass) return FALSE;
	SumMass = SumMom.m2();
	if (SumMass < 0) return FALSE;
	SumMass = std::sqrt(SumMass);

	// Compute c.m.s. hadron velocity and boost KTV to hadron c.m.s.
	G4ThreeVector Beta = -SumMom.boostVector();
	      //G4cout << " == pre boost 2 "<< SumMom.e()<< " "<< SumMom.mag()<<" "<< Beta <<G4endl;
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
			if (debug_G4BinaryLightIonReactionResults) G4cout << "E/p corrector: " << cAttempt << G4endl;
			success = true;
			break;
		}
		if ( cAttempt > 10 )
		{
			//         G4cout << " speed it up? " << std::abs(OldScale/(OldScale-Scale)) << G4endl;
			factor=std::max(1.,G4Log(std::abs(OldScale/(OldScale-Scale))));
			//	 G4cout << " ? factor ? " << factor << G4endl;
		}
	}

	if( (!success)  && debug_G4BinaryLightIonReactionResults)
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
G4bool G4BinaryLightIonReaction::SetLighterAsProjectile(G4LorentzVector & mom,const G4LorentzRotation & toBreit)
{
   G4bool swapped = false;
   if(tA<pA)
   {
      swapped = true;
      G4int tmp(0);
      tmp = tA; tA=pA; pA=tmp;
      tmp = tZ; tZ=pZ; pZ=tmp;
      G4double m1=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(pZ,pA);
      G4LorentzVector it(m1, G4ThreeVector(0,0,0));
      mom = toBreit*it;
   }
   return swapped;
}
G4ReactionProductVector * G4BinaryLightIonReaction::FuseNucleiAndPrompound(const G4LorentzVector & mom)
{
   // Check if kinematically nuclei can fuse.
   G4double mFused=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(pZ+tZ,pA+tA);
   G4double mTarget=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(tZ,tA);
   G4LorentzVector pCompound(mom.e()+mTarget,mom.vect());
   G4double m2Compound=pCompound.m2();
   if (m2Compound < sqr(mFused) ) {
	   //G4cout << "G4BLIC: projectile p, mTarget, mFused, mCompound, delta: " <<mom <<  " " << mTarget <<  " " << mFused
	   //	   <<  " " << sqrt(m2Compound)<<  " " << sqrt(m2Compound) - mFused << G4endl;
	   return 0;
   }

   G4Fragment aPreFrag;
   aPreFrag.SetZandA_asInt(pZ+tZ, pA+tA);
   aPreFrag.SetNumberOfParticles(pA);
   aPreFrag.SetNumberOfCharged(pZ);
   aPreFrag.SetNumberOfHoles(0);
   //GF FIXME: whyusing plop in z direction? this will not conserve momentum?
   //G4ThreeVector plop(0.,0., mom.vect().mag());
   //G4LorentzVector aL(mom.t()+mTarget, plop);
   G4LorentzVector aL(mom.t()+mTarget,mom.vect());
   aPreFrag.SetMomentum(aL);


         //G4cout << "Fragment INFO "<< pA+tA <<" "<<pZ+tZ<<" "
         //       << aL <<" "<<G4endl << aPreFrag << G4endl;
   G4ReactionProductVector * cascaders = theProjectileFragmentation->DeExcite(aPreFrag);
   //G4double tSum = 0;
   for(size_t count = 0; count<cascaders->size(); count++)
   {
      cascaders->operator[](count)->SetNewlyAdded(true);
      //tSum += cascaders->operator[](count)->GetKineticEnergy();
   }
   //       G4cout << "Exiting pre-compound only, E= "<<tSum<<G4endl;
   return cascaders;
}
G4ReactionProductVector * G4BinaryLightIonReaction::Interact(G4LorentzVector & mom, const G4LorentzRotation & toBreit)
{
      G4ReactionProductVector * result = 0;
      G4double projectileMass(0);
      G4LorentzVector it;

      G4int tryCount(0);
      do
      {
         ++tryCount;
         projectile3dNucleus = new G4Fancy3DNucleus;
         projectile3dNucleus->Init(pA, pZ);
         projectile3dNucleus->CenterNucleons();
         projectileMass=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(
               projectile3dNucleus->GetCharge(),projectile3dNucleus->GetMassNumber());
         it=toBreit * G4LorentzVector(projectileMass,G4ThreeVector(0,0,0));

         target3dNucleus = new G4Fancy3DNucleus;
         target3dNucleus->Init(tA, tZ);
         G4double impactMax = target3dNucleus->GetOuterRadius()+projectile3dNucleus->GetOuterRadius();
         //        G4cout << "out radius - nucleus - projectile " << target3dNucleus->GetOuterRadius()/fermi << " - " << projectile3dNucleus->GetOuterRadius()/fermi << G4endl;
         G4double aX=(2.*G4UniformRand()-1.)*impactMax;
         G4double aY=(2.*G4UniformRand()-1.)*impactMax;
         G4ThreeVector pos(aX, aY, -2.*impactMax-5.*fermi);

         G4KineticTrackVector * initalState = new G4KineticTrackVector;
         projectile3dNucleus->StartLoop();
         G4Nucleon * aNuc;
         G4LorentzVector tmpV(0,0,0,0);
         #ifdef debug_BLIR_finalstate
             G4LorentzVector pinitial;
         #endif
         G4LorentzVector nucleonMom(1./pA*mom);
         nucleonMom.setZ(nucleonMom.vect().mag());
         nucleonMom.setX(0);
         nucleonMom.setY(0);
         theFermi.Init(pA,pZ);
         while( (aNuc=projectile3dNucleus->GetNextNucleon()) )   /* Loop checking, 31.08.2015, G.Folger */
         {
            G4LorentzVector p4 = aNuc->GetMomentum();
            tmpV+=p4;
            G4ThreeVector nucleonPosition(aNuc->GetPosition());
            G4double density=(projectile3dNucleus->GetNuclearDensity())->GetDensity(nucleonPosition);
            nucleonPosition += pos;
            G4KineticTrack * it1 = new G4KineticTrack(aNuc, nucleonPosition, nucleonMom );
            it1->SetState(G4KineticTrack::outside);
            G4double pfermi= theFermi.GetFermiMomentum(density);
            G4double mass = aNuc->GetDefinition()->GetPDGMass();
            G4double Efermi= std::sqrt( sqr(mass) + sqr(pfermi)) - mass;
            it1->SetProjectilePotential(-Efermi);
            initalState->push_back(it1);
            #ifdef debug_BLIR_finalstate
               pinitial += it1->Get4Momentum();
            #endif
         }

         result=theModel->Propagate(initalState, target3dNucleus);
         #ifdef debug_BLIR_finalstate
         if( result && result->size()>0)
         {
        	 G4LorentzVector presult;
        	 G4ReactionProductVector::iterator iter;
        	 G4ReactionProduct xp;
        	 for (iter=result->begin(); iter !=result->end(); ++iter)
        	 {
        		presult += G4LorentzVector((*iter)->GetMomentum(),(*iter)->GetTotalEnergy());
        	 }

        	      G4cout << "BLIC check result :  initial " << pinitial << " mass tgt " << target3dNucleus->GetMass()
        			 << " final " << presult
        			 << " IF - FF " << pinitial +G4LorentzVector(target3dNucleus->GetMass()) - presult << G4endl;

         }
         #endif
         if( result && result->size()==0)
         {
            delete result;
            result=0;
         }
         if ( ! result )
         {
            delete target3dNucleus;
            delete projectile3dNucleus;
         }

         // std::for_each(initalState->begin(), initalState->end(), Delete<G4KineticTrack>());
         // delete initalState;

      } while (! result && tryCount< 150);   /* Loop checking, 31.08.2015, G.Folger */
      return result;
}
G4double G4BinaryLightIonReaction::GetProjectileExcitation()
{

      G4Nucleon * aNuc;
      // the projectileNucleus excitation energy estimate...
      G4double theStatisticalExEnergy = 0;
      projectile3dNucleus->StartLoop();
      while( (aNuc=projectile3dNucleus->GetNextNucleon()) )   /* Loop checking, 31.08.2015, G.Folger */
      {
                //G4cout << " Nucleon : " << aNuc->GetDefinition()->GetParticleName() <<" "<< aNuc->AreYouHit() <<" "<<aNuc->GetMomentum()<<G4endl;
         if(aNuc->AreYouHit()) {
            G4ThreeVector aPosition(aNuc->GetPosition());
            G4double localDensity = projectile3dNucleus->GetNuclearDensity()->GetDensity(aPosition);
            G4double localPfermi = theFermi.GetFermiMomentum(localDensity);
            G4double nucMass = aNuc->GetDefinition()->GetPDGMass();
            G4double localFermiEnergy = std::sqrt(nucMass*nucMass + localPfermi*localPfermi) - nucMass;
            G4double deltaE = localFermiEnergy - (aNuc->GetMomentum().t()-aNuc->GetMomentum().mag());
            theStatisticalExEnergy += deltaE;
         }
      }
      return theStatisticalExEnergy;
}

G4LorentzVector G4BinaryLightIonReaction::SortResult(G4ReactionProductVector * result, G4ReactionProductVector * spectators,G4ReactionProductVector * cascaders)
{
   unsigned int i(0);
   spectatorA=spectatorZ=0;
   G4LorentzVector pspectators(0,0,0,0);
   pFinalState=G4LorentzVector(0,0,0,0);
   for(i=0; i<result->size(); i++)
   {
      if( (*result)[i]->GetNewlyAdded() )
      {
         pFinalState += G4LorentzVector( (*result)[i]->GetMomentum(), (*result)[i]->GetTotalEnergy() );
         cascaders->push_back((*result)[i]);
      }
      else {
         //          G4cout <<" spectator ... ";
         pspectators += G4LorentzVector( (*result)[i]->GetMomentum(), (*result)[i]->GetTotalEnergy() );
         spectators->push_back((*result)[i]);
         spectatorA++;
         spectatorZ+= G4lrint((*result)[i]->GetDefinition()->GetPDGCharge()/eplus);
      }

      //       G4cout << (*result)[i]<< " "
      //        << (*result)[i]->GetDefinition()->GetParticleName() << " "
      //        << (*result)[i]->GetMomentum()<< " "
      //        << (*result)[i]->GetTotalEnergy() << G4endl;
   }
      //G4cout << "pFinalState / pspectators, (A,Z), p " << pFinalState << " / " << spectators->size()
      //		<< " (" << spectatorA << ", "<< spectatorZ  << "), 4-mom: " << pspectators << G4endl;

   return pspectators;
}

void G4BinaryLightIonReaction::DeExciteSpectatorNucleus(G4ReactionProductVector * spectators, G4ReactionProductVector * cascaders,
                                                 G4double theStatisticalExEnergy, G4LorentzVector & pSpectators)
{
   // call precompound model
   G4ReactionProductVector * proFrag = 0;
   G4LorentzVector pFragment(0.,0.,0.,0.);
   //      G4cout << " == pre boost 1 "<< momentum.e()<< " "<< momentum.mag()<<G4endl;
   G4LorentzRotation boost_fragments;
   //      G4cout << " == post boost 1 "<< momentum.e()<< " "<< momentum.mag()<<G4endl;
   //    G4LorentzRotation boost_spectator_mom(-momentum.boostVector());
   //     G4cout << "- momentum " << boost_spectator_mom * momentum << G4endl;
   G4LorentzVector pFragments(0,0,0,0);

   if(spectatorZ>0 && spectatorA>1)
   {
      //  Make the fragment
      G4Fragment aProRes;
      aProRes.SetZandA_asInt(spectatorZ, spectatorA);
      aProRes.SetNumberOfParticles(0);
      aProRes.SetNumberOfCharged(0);
      aProRes.SetNumberOfHoles(pA-spectatorA);
      G4double mFragment=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(spectatorZ,spectatorA);
      pFragment=G4LorentzVector(0,0,0,mFragment+std::max(0.,theStatisticalExEnergy) );
      aProRes.SetMomentum(pFragment);

      proFrag = theHandler->BreakItUp(aProRes);

      boost_fragments = G4LorentzRotation(pSpectators.boostVector());

      //     G4cout << " Fragment a,z, Mass Fragment, mass spect-mom, exitationE "
      //       << spectatorA <<" "<< spectatorZ <<" "<< mFragment <<" "
      //       << momentum.mag() <<" "<< momentum.mag() - mFragment
      //       << " "<<theStatisticalExEnergy
      //       << " "<< boost_fragments*pFragment<< G4endl;
      G4ReactionProductVector::iterator ispectator;
      for (ispectator=spectators->begin();ispectator!=spectators->end();ispectator++)
      {
         delete *ispectator;
      }
   }
   else if(spectatorA!=0)
   {
     G4ReactionProductVector::iterator ispectator;
     for (ispectator=spectators->begin();ispectator!=spectators->end();ispectator++)
      {
         (*ispectator)->SetNewlyAdded(true);
         cascaders->push_back(*ispectator);
         pFinalState+=G4LorentzVector((*ispectator)->GetMomentum(),(*ispectator)->GetTotalEnergy());
                  //G4cout << "BLIC: spectatorA>0, Z=0 from spectator "
                  // << (*ispectator)->GetDefinition()->GetParticleName() << " "
                  // << (*ispectator)->GetMomentum()<< " "
                  // << (*ispectator)->GetTotalEnergy() << G4endl;
      }

   }
   // / if (spectators)
   delete spectators;

   // collect the evaporation part and boost to spectator frame
   G4ReactionProductVector::iterator ii;
   if(proFrag)
   {
      for(ii=proFrag->begin(); ii!=proFrag->end(); ii++)
      {
         (*ii)->SetNewlyAdded(true);
         G4LorentzVector tmp((*ii)->GetMomentum(),(*ii)->GetTotalEnergy());
         tmp *= boost_fragments;
         (*ii)->SetMomentum(tmp.vect());
         (*ii)->SetTotalEnergy(tmp.e());
         //      result->push_back(*ii);
         pFragments += tmp;
      }
   }

   //    G4cout << "Fragmented p, momentum, delta " << pFragments <<" "<<momentum
   //            <<" "<< pFragments-momentum << G4endl;

   //  correct p/E of Cascade secondaries
   G4LorentzVector pCas=pInitialState - pFragments;

       //G4cout <<"BLIC: Going to correct from " << pFinalState << " to " << pCas << G4endl;
   //  the creation of excited fragment did violate E/p, so correct cascaders to get overall conservation.
   G4bool EnergyIsCorrect=EnergyAndMomentumCorrector(cascaders, pCas);
   if ( ! EnergyIsCorrect && debug_G4BinaryLightIonReactionResults)
   {
      G4cout << "G4BinaryLightIonReaction E/P correction for nucleus failed, will try to correct overall" << G4endl;
   }

   //  Add deexcitation secondaries
   if(proFrag)
   {
      for(ii=proFrag->begin(); ii!=proFrag->end(); ii++)
      {
         cascaders->push_back(*ii);
      }
      delete proFrag;
   }
      //G4cout << "EnergyIsCorrect? " << EnergyIsCorrect << G4endl;
   if ( ! EnergyIsCorrect )
   {
         // G4cout <<" ! EnergyIsCorrect " << pFinalState << " to " << pInitialState << G4endl;
      if (! EnergyAndMomentumCorrector(cascaders,pInitialState))
      {
         if(debug_G4BinaryLightIonReactionResults)
            G4cout << "G4BinaryLightIonReaction E/P corrections failed" << G4endl;
      }
   }

}

