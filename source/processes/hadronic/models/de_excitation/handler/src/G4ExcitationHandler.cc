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
// $Id: G4ExcitationHandler.cc 104780 2017-06-16 09:23:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
//
//
// Modified:
// 30 June 1998 by V. Lara:
//      -Modified the Transform method for use G4ParticleTable and 
//       therefore G4IonTable. It makes possible to convert all kind 
//       of fragments (G4Fragment) produced in deexcitation to 
//       G4DynamicParticle
//      -It uses default algorithms for:
//              Evaporation: G4Evaporation
//              MultiFragmentation: G4StatMF 
//              Fermi Breakup model: G4FermiBreakUp
// 24 Jul 2008 by M. A. Cortes Giraldo:
//      -Max Z,A for Fermi Break-Up turns to 9,17 by default
//      -BreakItUp() reorganised and bug in Evaporation loop fixed
//      -Transform() optimised
// (September 2008) by J. M. Quesada. External choices have been added for :
//      -inverse cross section option (default OPTxs=3)
//      -superimposed Coulomb barrier (if useSICB is set true, by default it is false) 
// September 2009 by J. M. Quesada: 
//      -according to Igor Pshenichnov, SMM will be applied (just in case) only once.
// 27 Nov 2009 by V.Ivanchenko: 
//      -cleanup the logic, reduce number internal vectors, fixed memory leak.
// 11 May 2010 by V.Ivanchenko: 
//      -FermiBreakUp activated, used integer Z and A, used BreakUpFragment method for 
//       final photon deexcitation; used check on adundance of a fragment, decay 
//       unstable fragments with A <5
// 22 March 2011 by V.Ivanchenko: general cleanup and addition of a condition: 
//       products of Fermi Break Up cannot be further deexcited by this model 
// 30 March 2011 by V.Ivanchenko removed private inline methods, moved Set methods 
//       to the source
// 23 January 2012 by V.Ivanchenko general cleanup including destruction of 
//    objects, propagate G4PhotonEvaporation pointer to G4Evaporation class and 
//    not delete it here 

#include "G4ExcitationHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4LorentzVector.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4Electron.hh"

#include "G4VMultiFragmentation.hh"
#include "G4VFermiBreakUp.hh"

#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Evaporation.hh"
#include "G4PhotonEvaporation.hh"
#include "G4StatMF.hh"
#include "G4FermiBreakUpVI.hh"
#include "G4NuclearLevelData.hh"
#include "G4Pow.hh"

G4ExcitationHandler::G4ExcitationHandler()
  : maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),
    fVerbose(0),isInitialised(false),isEvapLocal(true)
{                                                                          
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
  nist = G4NistManager::Instance();
  
  theMultiFragmentation = nullptr;
  theFermiModel = nullptr;
  G4Pow::GetInstance();
  theEvaporation = new G4Evaporation();
  thePhotonEvaporation = theEvaporation->GetPhotonEvaporation();
  theResults.reserve(60);
  results.reserve(30);
  theEvapList.reserve(30);
  thePhotoEvapList.reserve(10);

  SetParameters();
  electron = G4Electron::Electron();
  
  if(fVerbose > 0) { G4cout << "### New handler " << this << G4endl; }
}

G4ExcitationHandler::~G4ExcitationHandler()
{
  //G4cout << "### Delete handler " << this << G4endl;
  delete theMultiFragmentation;
  delete theFermiModel;
  if(isEvapLocal) { delete theEvaporation; } 
}

void G4ExcitationHandler::SetParameters()
{
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  minEForMultiFrag = param->GetMinExPerNucleounForMF();
  minExcitation = param->GetMinExcitation();
  icID = param->GetInternalConversionID();
  if(!thePhotonEvaporation) { 
    SetPhotonEvaporation(new G4PhotonEvaporation()); 
  } 
  if(!theFermiModel) { SetFermiModel(new G4FermiBreakUpVI()); }
  if(!theMultiFragmentation) { SetMultiFragmentation(new G4StatMF()); }
}

void G4ExcitationHandler::Initialise()
{
  if(isInitialised) { return; }
  if(fVerbose > 0) {
    G4cout << "G4ExcitationHandler::Initialise() started " << this << G4endl;
  }
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  isInitialised = true;
  SetParameters();
  theFermiModel->Initialise();
  theEvaporation->InitialiseChannels();
  param->Dump();
  if(G4Threading::IsMasterThread()) {
    G4cout << "Number of de-excitation channels                    " 
	   << theEvaporation->GetNumberOfChannels();
    if(fVerbose > 0) { G4cout << " " << this; }
    G4cout << G4endl; 
  }
}

void G4ExcitationHandler::SetEvaporation(G4VEvaporation* ptr, G4bool flag)
{
  if(ptr && ptr != theEvaporation) {
    delete theEvaporation; 
    theEvaporation = ptr;
    thePhotonEvaporation = ptr->GetPhotonEvaporation();
    theEvaporation->SetFermiBreakUp(theFermiModel);
    isEvapLocal = flag;
  }
}

void 
G4ExcitationHandler::SetMultiFragmentation(G4VMultiFragmentation* ptr)
{
  if(ptr && ptr != theMultiFragmentation) {
    delete theMultiFragmentation;
    theMultiFragmentation = ptr;
  }
}

void G4ExcitationHandler::SetFermiModel(G4VFermiBreakUp* ptr)
{
  if(ptr && ptr != theFermiModel) {
    delete theFermiModel;
    theFermiModel = ptr;
    theEvaporation->SetFermiBreakUp(theFermiModel);
  }
}

void 
G4ExcitationHandler::SetPhotonEvaporation(G4VEvaporationChannel* ptr)
{
  if(ptr && ptr != thePhotonEvaporation) {
    thePhotonEvaporation = ptr;
    theEvaporation->SetPhotonEvaporation(ptr);
  }
}

void G4ExcitationHandler::SetDeexChannelsType(G4DeexChannelType val)
{
  G4Evaporation* evap = static_cast<G4Evaporation*>(theEvaporation);
  if(!evap) { return; }
  if(val == fEvaporation) {
    evap->SetDefaultChannel();
  } else if(val == fCombined) {
    evap->SetCombinedChannel();
  } else if(val == fGEM) {
    evap->SetGEMChannel();
  }
  evap->InitialiseChannels();
  if(G4Threading::IsMasterThread()) {
    G4cout << "Number of de-excitation channels is changed to      " 
	   << theEvaporation->GetNumberOfChannels();
    if(fVerbose > 0) { G4cout << " " << this; }
    G4cout << G4endl; 
  }
}

G4ReactionProductVector * 
G4ExcitationHandler::BreakItUp(const G4Fragment & theInitialState)
{
  // Variables existing until end of method
  G4Fragment * theInitialStatePtr = new G4Fragment(theInitialState);
  if(fVerbose > 1) { 	
    G4cout << "@@@@@@@@@@ Start G4Excitation Handler @@@@@@@@@@@@@ " << G4endl;
    G4cout << theInitialState << G4endl;  
  }
  if(!isInitialised) { Initialise(); }

  // pointer to fragment vector which receives temporal results
  G4FragmentVector * theTempResult = nullptr;

  theResults.clear();
  thePhotoEvapList.clear();
  theEvapList.clear();
   
  // Variables to describe the excited configuration
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4int A = theInitialState.GetA_asInt();
  G4int Z = theInitialState.GetZ_asInt();
  
  // In case A <= 1 the fragment will not perform any nucleon emission
  if (A <= 1) {
    theResults.push_back( theInitialStatePtr );

    // check if a fragment is stable
  } else if(exEnergy < minExcitation && nist->GetIsotopeAbundance(Z, A) > 0.0) {
    theResults.push_back( theInitialStatePtr );

    // JMQ 150909: first step in de-excitation is treated separately 
    // Fragments after the first step are stored in theEvapList 
  } else {      
    if((A<maxAForFermiBreakUp && Z<maxZForFermiBreakUp) 
       || exEnergy <= minEForMultiFrag*A) { 
      theEvapList.push_back(theInitialStatePtr); 

    // Statistical Multifragmentation will take place only once
    } else {
      theTempResult = theMultiFragmentation->BreakItUp(theInitialState);
      if(!theTempResult) { 
	theEvapList.push_back(theInitialStatePtr); 
      } else {
	size_t nsec = theTempResult->size();

	// no fragmentation
	if(0 == nsec) { 
	  theEvapList.push_back(theInitialStatePtr); 

	  // secondary are produced - sort out secondary fragments
	} else {
	  G4bool deletePrimary = true;
	  G4FragmentVector::iterator j;
	  for (j = theTempResult->begin(); j != theTempResult->end(); ++j) {  

	    if((*j) == theInitialStatePtr) { deletePrimary = false; }
	    A = (*j)->GetA_asInt();  

	    // gamma, p, n
	    if(A <= 1) { 
	      theResults.push_back(*j); 

	      // Analyse fragment A > 1
	    } else {
	      G4double exEnergy1 = (*j)->GetExcitationEnergy();

	      // cold fragments
	      if(exEnergy1 < minExcitation) {
		Z = (*j)->GetZ_asInt(); 
		if(nist->GetIsotopeAbundance(Z, A) > 0.0) { 
		  theResults.push_back(*j); // stable fragment 
		} else {
		  theEvapList.push_back(*j);
		}
		// hot fragments are unstable
	      } else { 
		theEvapList.push_back(*j); 
	      } 
	    }
	  }
	  if( deletePrimary ) { delete theInitialStatePtr; }
	}
	delete theTempResult; // end multifragmentation
      }
    }
  }
  if(fVerbose > 2) { 	
    G4cout << "## After first step " << theEvapList.size() << " for evap;  "
	   << thePhotoEvapList.size() << " for photo-evap; " 
	   << theResults.size() << " results. " << G4endl; 
  }
  // -----------------------------------
  // FermiBreakUp and De-excitation loop
  // -----------------------------------
      
  static const G4int countmax = 1000;
  G4Fragment* frag;
  size_t kk;
  for (kk=0; kk<theEvapList.size(); ++kk) {
    frag = theEvapList[kk];
    if(fVerbose > 2) { 	
      G4cout << "Next evaporate: " << G4endl;  
      G4cout << *frag << G4endl;  
    }
    if(kk >= countmax) {
      G4ExceptionDescription ed;
      ed << "Infinite loop in the de-excitation module: " << kk
	 << " iterations \n"
	 << "      Initial fragment: \n" << theInitialState
	 << "\n      Current fragment: \n" << *frag;
      G4Exception("G4ExcitationHandler::BreakItUp","had0333",FatalException,
		  ed,"Stop execution");
      
    }
    A = frag->GetA_asInt(); 
    Z = frag->GetZ_asInt();
    results.clear();
	  
    // Fermi Break-Up 
    if(theFermiModel->IsApplicable(Z, A, frag->GetExcitationEnergy())) {
      theFermiModel->BreakFragment(&results, frag);
      size_t nsec = results.size();
      if(fVerbose > 2) { G4cout << "FermiBreakUp Nsec= " << nsec << G4endl; }

      // FBU takes care to delete input fragment or add it to the results
      // The secondary may be excited - photo-evaporation should be applied
      for(size_t j=0; j<nsec; ++j) {
	exEnergy = results[j]->GetExcitationEnergy();
	if(exEnergy < minExcitation) { theResults.push_back(results[j]); }
	else                   { thePhotoEvapList.push_back(results[j]); }
      }
      continue; 
    }
    // apply Evaporation, residual nucleus is always added to the results
    theEvaporation->BreakFragment(&results, frag); 
    size_t nsec = results.size();
    if(fVerbose > 2) { G4cout << "Evaporation Nsec= " << nsec << G4endl; }
		  
    // no evaporation
    if(1 >= nsec) {
      theResults.push_back(frag);
      continue;
    }

    // Sort out secondary fragments
    for (size_t j = 0; j<nsec; ++j) {
      if(fVerbose > 3) {
	G4cout << "Evaporated product #" << j << G4endl;
	G4cout << results[j] << G4endl;  
      }
      A = results[j]->GetA_asInt();
      if(A <= 1) { 
	theResults.push_back(results[j]);   // gamma, p, n
	continue;
      } 
      exEnergy = results[j]->GetExcitationEnergy();

      // hot fragment
      if(exEnergy >= minExcitation) {
	theEvapList.push_back(results[j]);        
	    
	// cold fragment
      } else {  
	Z = results[j]->GetZ_asInt();

	// natural isotope
	if(nist->GetIsotopeAbundance(Z, A) > 0.0) { 
	  theResults.push_back(results[j]); // stable fragment 

	} else {
	  theEvapList.push_back(results[j]);   	  
	}
      }
    } // end of loop on secondary
  } // end of the loop over theEvapList
  if(fVerbose > 2) { 	
    G4cout << "## After 2nd step " << theEvapList.size() << " was evap;  "
	   << thePhotoEvapList.size() << " for photo-evap; " 
	   << theResults.size() << " results. " << G4endl; 
  }
  // -----------------------
  // Photon-Evaporation loop
  // -----------------------
  
  // at this point only photon evaporation is possible
  size_t kkmax = thePhotoEvapList.size();
  for (kk=0; kk<kkmax; ++kk) {
    frag = thePhotoEvapList[kk];
    if(fVerbose > 2) { 	
      G4cout << "Next photon evaporate: " << thePhotonEvaporation << G4endl;  
      G4cout << *frag << G4endl;
    }
    exEnergy = frag->GetExcitationEnergy();

    // photon de-excitation only for hot fragments
    if(exEnergy > minExcitation) {  
      thePhotonEvaporation->BreakUpChain(&theResults, frag);
    }

    // primary fragment is kept
    theResults.push_back(frag); 

  } // end of photon-evaporation loop
  if(fVerbose > 2) { 	
    G4cout << "## After 3d step " << theEvapList.size() << " was evap;  "
	   << thePhotoEvapList.size() << " was photo-evap; " 
	   << theResults.size() << " results. " << G4endl; 
  }
  G4ReactionProductVector * theReactionProductVector = 
    new G4ReactionProductVector();

  // MAC (24/07/08)
  // To optimise the storing speed, we reserve space in memory for the vector
  theReactionProductVector->reserve( theResults.size() );

  G4int theFragmentA, theFragmentZ;

  if(fVerbose > 1) { 	
    G4cout << "### ExcitationHandler provides " << theResults.size() 
	   << " evaporated products:" << G4endl;
  }
  kkmax = theResults.size();
  for (kk=0; kk<kkmax; ++kk) {
    frag = theResults[kk];
    if(fVerbose > 1) { 
      G4cout << kk << "-th fragment " << frag;
      if(frag->GetNuclearPolarization()) { G4cout << "  " << frag->GetNuclearPolarization(); }
      G4cout << G4endl;
      G4cout << *frag << G4endl; 
    }

    theFragmentA = frag->GetA_asInt();
    theFragmentZ = frag->GetZ_asInt();
    G4double etot= frag->GetMomentum().e();
    G4double eexc = 0.0;
    const G4ParticleDefinition* theKindOfFragment = nullptr;
    if (theFragmentA == 0) {       // photon or e-
      theKindOfFragment = frag->GetParticleDefinition();   
    } else if (theFragmentA == 1 && theFragmentZ == 0) { // neutron
      theKindOfFragment = G4Neutron::NeutronDefinition();
    } else if (theFragmentA == 1 && theFragmentZ == 1) { // proton
      theKindOfFragment = G4Proton::ProtonDefinition();
    } else if (theFragmentA == 2 && theFragmentZ == 1) { // deuteron
      theKindOfFragment = G4Deuteron::DeuteronDefinition();
    } else if (theFragmentA == 3 && theFragmentZ == 1) { // triton
      theKindOfFragment = G4Triton::TritonDefinition();
    } else if (theFragmentA == 3 && theFragmentZ == 2) { // helium3
      theKindOfFragment = G4He3::He3Definition();
    } else if (theFragmentA == 4 && theFragmentZ == 2) { // alpha
      theKindOfFragment = G4Alpha::AlphaDefinition();;
    } else {

      // fragment
      eexc = frag->GetExcitationEnergy();
      G4int idxf = frag->GetFloatingLevelNumber();
      if(eexc < minExcitation) { 
	eexc = 0.0; 
        idxf = 0;
      }

      theKindOfFragment = theTableOfIons->GetIon(theFragmentZ,theFragmentA,eexc,
                                                 G4Ions::FloatLevelBase(idxf));
      if(fVerbose > 1) {
	G4cout << "### EXCH: Find ion Z= " << theFragmentZ << " A= " << theFragmentA
	       << " Eexc(MeV)= " << eexc/MeV << " idx= " << idxf 
	       << "  " <<  theKindOfFragment << G4endl;
      }
    }
    // fragment identified
    if(theKindOfFragment) {
      G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
      theNew->SetMomentum(frag->GetMomentum().vect());
      theNew->SetTotalEnergy(etot);
      theNew->SetFormationTime(frag->GetCreationTime());
      if(theKindOfFragment == electron) { theNew->SetCreatorModel(icID); }
      theReactionProductVector->push_back(theNew);

      // fragment not found out ground state is created
    } else { 
      theKindOfFragment = 
	theTableOfIons->GetIon(theFragmentZ,theFragmentA,0.0,noFloat,0);
      if(theKindOfFragment) {
	G4ThreeVector mom(0.0,0.0,0.0); 
	G4double ionmass = theKindOfFragment->GetPDGMass();
	if(etot <= ionmass) {
	  etot = ionmass;
	} else {
	  G4double ptot = std::sqrt((etot - ionmass)*(etot + ionmass));
	  mom = (frag->GetMomentum().vect().unit())*ptot;
	}
	G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
	theNew->SetMomentum(mom);
	theNew->SetTotalEnergy(etot);
	theNew->SetFormationTime(frag->GetCreationTime());
	theReactionProductVector->push_back(theNew);
	if(fVerbose > 2) {            	
	  G4cout << "### Find ion Z= " << theFragmentZ << " A= " << theFragmentA
		 << " ground state, energy corrected E(MeV)= " << etot << G4endl;
	}
      }
    }
    delete frag;
    if(fVerbose > 1) { G4cout << "G4Fragment #" << kk << " is deleted" << G4endl; }
  }
  if(fVerbose > 2) { 	
    G4cout << "@@@@@@@@@@ End G4Excitation Handler "<< G4endl;
  }
  return theReactionProductVector;
}

void G4ExcitationHandler::ModelDescription(std::ostream& outFile) const
{
  outFile << "G4ExcitationHandler description\n"
	  << "This class samples de-excitation of excited nucleus using\n"
	  << "Fermi Break-up model for light fragments (Z < 9, A < 17), "
	  << "evaporation, fission, and photo-evaporation models. Evaporated\n"
	  << "particle may be proton, neutron, and other light fragment \n"
	  << "(Z < 13, A < 29). During photon evaporation produced gamma \n"
	  << "or electrons due to internal conversion \n";
}






