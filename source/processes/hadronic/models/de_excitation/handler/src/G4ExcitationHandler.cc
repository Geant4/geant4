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
#include "G4ThreeVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4Electron.hh"
#include "G4Lambda.hh"

#include "G4VMultiFragmentation.hh"
#include "G4VFermiBreakUp.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4Evaporation.hh"
#include "G4PhotonEvaporation.hh"
#include "G4StatMF.hh"
#include "G4FermiBreakUpVI.hh"
#include "G4NuclearLevelData.hh"
#include "G4Pow.hh"
#include "G4PhysicsModelCatalog.hh"

G4ExcitationHandler::G4ExcitationHandler()
  : icID(0),maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),
    fVerbose(1),fWarnings(0),minEForMultiFrag(1.*CLHEP::TeV),
    minExcitation(1.*CLHEP::eV),maxExcitation(100.*CLHEP::MeV),
    isInitialised(false),isEvapLocal(true),isActive(true)
{                                                                          
  thePartTable = G4ParticleTable::GetParticleTable();
  theTableOfIons = thePartTable->GetIonTable();
  nist = G4NistManager::Instance();
  
  theMultiFragmentation = nullptr;
  theFermiModel = nullptr;
  theEvaporation = nullptr;
  thePhotonEvaporation = nullptr;
  theResults.reserve(60);
  results.reserve(30);
  theEvapList.reserve(30);

  G4Pow::GetInstance();
  theElectron = G4Electron::Electron();
  theNeutron = G4Neutron::NeutronDefinition();
  theProton = G4Proton::ProtonDefinition();
  theDeuteron = G4Deuteron::DeuteronDefinition();
  theTriton = G4Triton::TritonDefinition();
  theHe3 = G4He3::He3Definition();
  theAlpha = G4Alpha::AlphaDefinition();
  theLambda = G4Lambda::Lambda();

  fLambdaMass = theLambda->GetPDGMass();

  if(fVerbose > 1) { G4cout << "### New handler " << this << G4endl; }
}

G4ExcitationHandler::~G4ExcitationHandler()
{
  delete theMultiFragmentation;
  delete theFermiModel;
  if(isEvapLocal) { delete theEvaporation; } 
}

void G4ExcitationHandler::SetParameters()
{
  G4NuclearLevelData* ndata = G4NuclearLevelData::GetInstance();
  auto param = ndata->GetParameters();
  isActive = true;
  // check if de-excitation is needed
  if(fDummy == param->GetDeexChannelsType()) { 
    isActive = false; 
  } else {
    // upload data for elements used in geometry
    G4int Zmax = 20;
    const G4ElementTable* table = G4Element::GetElementTable();
    for(auto & elm : *table) { Zmax = std::max(Zmax, elm->GetZasInt()); }
    ndata->UploadNuclearLevelData(Zmax+1);
  }
  minEForMultiFrag = param->GetMinExPerNucleounForMF();
  minExcitation = param->GetMinExcitation();
  maxExcitation = param->GetPrecoHighEnergy();
  icID = G4PhysicsModelCatalog::GetModelID("model_e-InternalConversion");

  // allowing local debug printout 
  fVerbose = std::max(fVerbose, param->GetVerbose());
  if(isActive) {
    if(!thePhotonEvaporation)  { SetPhotonEvaporation(new G4PhotonEvaporation()); }
    if(!theEvaporation) { 
      SetEvaporation(new G4Evaporation(thePhotonEvaporation), true); 
    }
    if(!theFermiModel)         { SetFermiModel(new G4FermiBreakUpVI()); }
    if(!theMultiFragmentation) { SetMultiFragmentation(new G4StatMF()); }
  }
  theFermiModel->SetVerbose(fVerbose);
  if(fVerbose > 1) {
    G4cout << "G4ExcitationHandler::SetParameters() done " << this << G4endl;
  }
}

void G4ExcitationHandler::Initialise()
{
  if(isInitialised) { return; }
  if(fVerbose > 1) {
    G4cout << "G4ExcitationHandler::Initialise() started " << this << G4endl;
  }
  G4DeexPrecoParameters* param = 
    G4NuclearLevelData::GetInstance()->GetParameters();
  isInitialised = true;
  SetParameters();
  if(isActive) {
    theFermiModel->Initialise();
    theEvaporation->InitialiseChannels();
  }
  // dump level is controlled by parameter class
  param->Dump();
}

void G4ExcitationHandler::SetEvaporation(G4VEvaporation* ptr, G4bool flag)
{
  if(ptr && ptr != theEvaporation) {
    delete theEvaporation; 
    theEvaporation = ptr;
    SetPhotonEvaporation(ptr->GetPhotonEvaporation());
    theEvaporation->SetFermiBreakUp(theFermiModel);
    isEvapLocal = flag;
    if(fVerbose > 1) {
      G4cout << "G4ExcitationHandler::SetEvaporation() for " << this << G4endl;
    }
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
    if(theEvaporation) { theEvaporation->SetFermiBreakUp(theFermiModel); }
  }
}

void 
G4ExcitationHandler::SetPhotonEvaporation(G4VEvaporationChannel* ptr)
{
  if(ptr && ptr != thePhotonEvaporation) {
    delete thePhotonEvaporation;
    thePhotonEvaporation = ptr;
    if(theEvaporation) { theEvaporation->SetPhotonEvaporation(ptr); }
    if(fVerbose > 1) {
      G4cout << "G4ExcitationHandler::SetPhotonEvaporation() " << ptr 
             << " for handler " << this << G4endl;
    }
  }
}

void G4ExcitationHandler::SetDeexChannelsType(G4DeexChannelType val)
{
  G4Evaporation* evap = static_cast<G4Evaporation*>(theEvaporation);
  if(fVerbose > 1) {
    G4cout << "G4ExcitationHandler::SetDeexChannelsType " << val 
	   << " for " << this << G4endl;
  }
  if(val == fDummy) { 
    isActive = false;
    return; 
  }
  if(!evap) { return; }
  if(val == fEvaporation) {
    evap->SetDefaultChannel();
  } else if(val == fCombined) {
    evap->SetCombinedChannel();
  } else if(val == fGEM) {
    evap->SetGEMChannel();
  } else if(val == fGEMVI) {
    evap->SetGEMVIChannel();
  }
  evap->InitialiseChannels();
  if(fVerbose > 1) {
    if(G4Threading::IsMasterThread()) {
      G4cout << "Number of de-excitation channels is changed to: " 
	     << theEvaporation->GetNumberOfChannels();
      G4cout << " " << this;
    }
    G4cout << G4endl; 
  }
}

G4VEvaporation* G4ExcitationHandler::GetEvaporation()
{
  if(!theEvaporation) { SetParameters(); }
  return theEvaporation;
}

G4VMultiFragmentation* G4ExcitationHandler::GetMultiFragmentation()
{
  if(!theMultiFragmentation) { SetParameters(); }
  return theMultiFragmentation;
}

G4VFermiBreakUp* G4ExcitationHandler::GetFermiModel()
{
  if(!theFermiModel) { SetParameters(); }
  return theFermiModel;
}

G4VEvaporationChannel* G4ExcitationHandler::GetPhotonEvaporation()
{
  if(!thePhotonEvaporation) { SetParameters(); }
  return thePhotonEvaporation;
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
  theEvapList.clear();
   
  // Variables to describe the excited configuration
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4int A = theInitialState.GetA_asInt();
  G4int Z = theInitialState.GetZ_asInt();
  G4int nL = theInitialState.GetNumberOfLambdas();

  // too much excitation
  if(exEnergy > A*maxExcitation && A > 0) {
    ++fWarnings;
    if(fWarnings < 0) {
      G4ExceptionDescription ed;
      ed << "High excitation Fragment Z= " << Z << " A= " << A 
	 << " Eex/A(MeV)= " << exEnergy/A;
      G4Exception("G4ExcitationHandler::BreakItUp()","had0034",JustWarning,ed,"");
    }
  }

  // for hyper-nuclei subtract lambdas from the projectile fragment
  G4double lambdaF = 0.0;
  G4LorentzVector lambdaLV = theInitialStatePtr->GetMomentum();
  if(0 < nL) {

    // is it a stable hyper-nuclei?
    if(A >= 3 && A <= 5 && nL <= 2) {
      G4int pdg = 0;
      if(3 == A && 1 == nL) {
        pdg = 1010010030;
      } else if(5 == A && 2 == Z && 1 == nL) {
        pdg = 1010020050;
      } else if(4 == A) {
	if(1 == Z && 1 == nL) {
	  pdg = 1010010040;
	} else if(2 == Z && 1 == nL) {
	  pdg = 1010020040;
	} else if(0 == Z && 2 == nL) {
	  pdg = 1020000040;
	} else if(1 == Z && 2 == nL) {
	  pdg = 1020010040;
	}
      }
      // initial state is one of hyper-nuclei
      if(0 < pdg) {
	const G4ParticleDefinition* part = thePartTable->FindParticle(pdg);
        if(nullptr != part) {
	  G4ReactionProduct* theNew = new G4ReactionProduct(part);
	  G4ThreeVector dir = G4ThreeVector( 0.0, 0.0, 0.0 );
          if ( lambdaLV.vect().mag() > CLHEP::eV ) {
	    dir = lambdaLV.vect().unit();
          }
	  G4double mass = part->GetPDGMass(); 
	  G4double etot = std::max(lambdaLV.e(), mass);
          dir *= std::sqrt((etot - mass)*(etot + mass));
	  theNew->SetMomentum(dir);
	  theNew->SetTotalEnergy(etot);
	  theNew->SetFormationTime(theInitialState.GetCreationTime());
	  theNew->SetCreatorModelID(theInitialState.GetCreatorModelID());
	  G4ReactionProductVector* v = new G4ReactionProductVector();
          v->push_back(theNew);
	  return v;
	}
      }
    }
    G4double mass = theInitialStatePtr->GetGroundStateMass();
    lambdaF = nL*(fLambdaMass - CLHEP::neutron_mass_c2)/mass;

    // de-excitation with neutrons instead of lambda inside the fragment
    theInitialStatePtr->SetZAandMomentum(lambdaLV*(1. - lambdaF), Z, A, 0);

    // 4-momentum not used in de-excitation 
    lambdaLV *= lambdaF;
  } else if(0 > nL) {
    ++fWarnings;
    if(fWarnings < 0) {
      G4ExceptionDescription ed;
      ed << "Fragment with negative L: Z=" << Z << " A=" << A << " L=" << nL
	 << " Eex/A(MeV)= " << exEnergy/A;
      G4Exception("G4ExcitationHandler::BreakItUp()","had0034",JustWarning,ed,"");
    }
  }

  // In case A <= 1 the fragment will not perform any nucleon emission
  if (A <= 1 || !isActive) {
    theResults.push_back( theInitialStatePtr );

    // check if a fragment is stable
  } else if(exEnergy < minExcitation && 
	    nist->GetIsotopeAbundance(Z, A) > 0.0) {
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
	  for (auto ptr : *theTempResult) {  
	    if(ptr == theInitialStatePtr) { deletePrimary = false; }
	    SortSecondaryFragment(ptr);
	  }
	  if( deletePrimary ) { delete theInitialStatePtr; }
	}
	delete theTempResult; // end multifragmentation
      }
    }
  }
  if(fVerbose > 2) { 	
    G4cout << "## After first step of handler " << theEvapList.size() 
           << " for evap;  "
	   << theResults.size() << " results. " << G4endl; 
  }
  // -----------------------------------
  // FermiBreakUp and De-excitation loop
  // -----------------------------------
      
  static const G4int countmax = 1000;
  size_t kk;
  for (kk=0; kk<theEvapList.size(); ++kk) {
    G4Fragment* frag = theEvapList[kk];
    if(fVerbose > 3) { 	
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
    if(fVerbose > 2) {
      G4cout << "G4ExcitationHandler# " << kk << " Z= " << Z << " A= " << A 
             << " Eex(MeV)= " << frag->GetExcitationEnergy() << G4endl;
    }	  
    // Fermi Break-Up 
    if(theFermiModel->IsApplicable(Z, A, frag->GetExcitationEnergy())) {
      theFermiModel->BreakFragment(&results, frag);
      size_t nsec = results.size();
      if(fVerbose > 2) { G4cout << "FermiBreakUp Nsec= " << nsec << G4endl; }

      // FBU takes care to delete input fragment or add it to the results
      // The secondary may be excited - photo-evaporation should be applied
      if(1 < nsec) {
	for(auto & res : results) {
	  SortSecondaryFragment(res);
	}
	continue;
      }
      // evaporation will be applied
    }
    // apply Evaporation, residual nucleus is always added to the results
    // photon evaporation is possible 
    theEvaporation->BreakFragment(&results, frag); 
    if(fVerbose > 3) { 
      G4cout << "Evaporation Nsec= " << results.size() << G4endl; 
    }
    if(0 == results.size()) {
      theResults.push_back(frag);
    } else {
      SortSecondaryFragment(frag);
    }

    // Sort out secondary fragments
    for (auto & res : results) {
      if(fVerbose > 4) {
	G4cout << "Evaporated product #" << *res << G4endl;  
      }
      SortSecondaryFragment(res);
    } // end of loop on secondary
  } // end of the loop over theEvapList
  if(fVerbose > 2) { 	
    G4cout << "## After 2nd step of handler " << theEvapList.size() 
           << " was evap;  "
	   << theResults.size() << " results. " << G4endl; 
  }
  G4ReactionProductVector * theReactionProductVector = 
    new G4ReactionProductVector();

  // MAC (24/07/08)
  // To optimise the storing speed, we reserve space 
  // in memory for the vector
  theReactionProductVector->reserve( theResults.size() );

  if(fVerbose > 2) { 	
    G4cout << "### ExcitationHandler provides " << theResults.size() 
	   << " evaporated products:" << G4endl;
  }
  G4LorentzVector partOfLambdaLV;
  if ( nL > 0 ) partOfLambdaLV = lambdaLV/(G4double)nL;
  for (auto & frag : theResults) {
    G4LorentzVector lv0 = frag->GetMomentum();
    G4double etot = lv0.e();

    // in the case of dummy de-excitation, excitation energy is transfered 
    // into kinetic energy of output ion
    if(!isActive) {
      G4double mass = frag->GetGroundStateMass();
      G4double ptot = lv0.vect().mag();
      G4double fac  = (etot <= mass || 0.0 == ptot) ? 0.0 
	: std::sqrt((etot - mass)*(etot + mass))/ptot; 
      G4LorentzVector lv((frag->GetMomentum()).px()*fac, 
			 (frag->GetMomentum()).py()*fac,
			 (frag->GetMomentum()).pz()*fac, etot);
      frag->SetMomentum(lv);
    }
    if(fVerbose > 3) { 
      G4cout << *frag;
      if(frag->NuclearPolarization()) { 
	G4cout << "  " << frag->NuclearPolarization(); 
      }
      G4cout << G4endl;
    }

    G4int fragmentA = frag->GetA_asInt();
    G4int fragmentZ = frag->GetZ_asInt();
    G4double eexc = 0.0;
    const G4ParticleDefinition* theKindOfFragment = nullptr;
    G4bool isHyperN = false;
    if (fragmentA == 0) {       // photon or e-
      theKindOfFragment = frag->GetParticleDefinition();   
    } else if (fragmentA == 1 && fragmentZ == 0) { // neutron
      theKindOfFragment = theNeutron;
    } else if (fragmentA == 1 && fragmentZ == 1) { // proton
      theKindOfFragment = theProton;
    } else if (fragmentA == 2 && fragmentZ == 1) { // deuteron
      theKindOfFragment = theDeuteron;
    } else if (fragmentA == 3 && fragmentZ == 1) { // triton
      theKindOfFragment = theTriton;
      if(0 < nL) {
        const G4ParticleDefinition* p = thePartTable->FindParticle(1010010030);
        if(nullptr != p) {
	  theKindOfFragment = p;
	  isHyperN = true;
	  --nL;
	}
      }
    } else if (fragmentA == 3 && fragmentZ == 2) { // helium3
      theKindOfFragment = theHe3;
    } else if (fragmentA == 4 && fragmentZ == 2) { // alpha
      theKindOfFragment = theAlpha;
      if(0 < nL) {
        const G4ParticleDefinition* p = thePartTable->FindParticle(1010020040);
        if(nullptr != p) {
	  theKindOfFragment = p;
	  isHyperN = true;
	  --nL;
	}
      }
    } else {

      // fragment
      eexc = frag->GetExcitationEnergy();
      G4int idxf = frag->GetFloatingLevelNumber();
      if(eexc < minExcitation) {
	eexc = 0.0; 
        idxf = 0;
      }

      theKindOfFragment = theTableOfIons->GetIon(fragmentZ, fragmentA, eexc,
                                                 G4Ions::FloatLevelBase(idxf));
      if(fVerbose > 3) {
	G4cout << "### EXCH: Find ion Z= " << fragmentZ 
               << " A= " << fragmentA
	       << " Eexc(MeV)= " << eexc/MeV << " idx= " << idxf 
	       << G4endl;
      }
    }
    // fragment identified
    if(nullptr != theKindOfFragment) {
      G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
      if(isHyperN) {
        G4LorentzVector lv = lv0 + partOfLambdaLV;
	G4ThreeVector dir = lv.vect().unit();
        G4double mass = theKindOfFragment->GetPDGMass();
        etot = std::max(lv.e(), mass);
        G4double ptot = std::sqrt((etot - mass)*(etot + mass));
        dir *= ptot;
	theNew->SetMomentum(dir);
	// remaining not compensated 4-momentum
        lambdaLV += (lv0 - G4LorentzVector(dir, etot));
      } else {
	theNew->SetMomentum(lv0.vect());
      }
      theNew->SetTotalEnergy(etot);
      theNew->SetFormationTime(frag->GetCreationTime());
      if(theKindOfFragment == theElectron) {
	theNew->SetCreatorModelID(icID);
      } else {
	theNew->SetCreatorModelID(frag->GetCreatorModelID());
      }
      theReactionProductVector->push_back(theNew);

      // fragment not found out ground state is created
    } else { 
      theKindOfFragment = 
	theTableOfIons->GetIon(fragmentZ,fragmentA,0.0,noFloat,0);
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
	theNew->SetCreatorModelID(frag->GetCreatorModelID());
	theReactionProductVector->push_back(theNew);
	if(fVerbose > 3) {
	  G4cout << "          ground state, energy corrected E(MeV)= " 
                 << etot << G4endl;
	}
      }
    }
    delete frag;
  }
  // remaining lambdas are free; conserve quantum numbers but
  // not 4-momentum
  if(0 < nL) {
    G4ThreeVector dir = G4ThreeVector( 0.0, 0.0, 0.0 );
    if ( lambdaLV.vect().mag() > CLHEP::eV ) {
      dir = lambdaLV.vect().unit();
    }
    G4double etot = std::max(lambdaLV.e()/(G4double)nL, fLambdaMass);
    dir *= std::sqrt((etot - fLambdaMass)*(etot + fLambdaMass));
    for(G4int i=0; i<nL; ++i) {
      G4ReactionProduct* theNew = new G4ReactionProduct(theLambda);
      theNew->SetMomentum(dir);
      theNew->SetTotalEnergy(etot);
      theNew->SetFormationTime(theInitialState.GetCreationTime());
      theNew->SetCreatorModelID(theInitialState.GetCreatorModelID());
      theReactionProductVector->push_back(theNew);
    }
  }
  if(fVerbose > 3) {
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



