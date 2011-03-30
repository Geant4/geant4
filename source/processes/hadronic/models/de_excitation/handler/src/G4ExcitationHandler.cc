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
// $Id: G4ExcitationHandler.cc,v 1.40 2010-11-17 16:20:38 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
//
//
// Modified:
// (30 June 1998) by V. Lara:
//      -Modified the Transform method for use G4ParticleTable and 
//       therefore G4IonTable. It makes possible to convert all kind 
//       of fragments (G4Fragment) produced in deexcitation to 
//       G4DynamicParticle
//      -It uses default algorithms for:
//              Evaporation: G4Evaporation
//              MultiFragmentation: G4StatMF 
//              Fermi Breakup model: G4FermiBreakUp
// (24 Jul 2008) by M. A. Cortes Giraldo:
//      -Max Z,A for Fermi Break-Up turns to 9,17 by default
//      -BreakItUp() reorganised and bug in Evaporation loop fixed
//      -Transform() optimised
// (September 2008) by J. M. Quesada. External choices have been added for :
//      -inverse cross section option (default OPTxs=3)
//      -superimposed Coulomb barrier (if useSICB is set true, by default it is false) 
// (September 2009) by J. M. Quesada: 
//      -according to Igor Pshenichnov, SMM will be applied (just in case) only once.
// (27 Nov 2009) by V.Ivanchenko: 
//      -cleanup the logic, reduce number internal vectors, fixed memory leak.
// (11 May 2010) by V.Ivanchenko: 
//      -FermiBreakUp activated, used integer Z and A, used BreakUpFragment method for 
//       final photon deexcitation; used check on adundance of a fragment, decay 
//       unstable fragments with A <5
// (22 March 2011) by V.Ivanchenko: general cleanup and addition of a condition: 
//       products of Fermi Break Up cannot be further deexcited by this model 
// (30 March 2011) by V.Ivanchenko removed private inline methods, moved Set methods 
//       to the source

#include "G4ExcitationHandler.hh"
#include "G4LorentzVector.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"

#include "G4VMultiFragmentation.hh"
#include "G4VFermiBreakUp.hh"
#include "G4VEvaporation.hh"
#include "G4VEvaporationChannel.hh"
#include "G4VPhotonEvaporation.hh"
#include "G4Evaporation.hh"
#include "G4StatMF.hh"
#include "G4FermiBreakUp.hh"
#include "G4PhotonEvaporation.hh"

#include <list>

G4ExcitationHandler::G4ExcitationHandler():
  maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),minEForMultiFrag(4.0*GeV),
  //  maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),minEForMultiFrag(3.0*MeV),
  minExcitation(CLHEP::keV),MyOwnEvaporationClass(true),
  MyOwnPhotonEvaporationClass(true),OPTxs(3),useSICB(false)
{                                                                          
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
  
  theEvaporation = new G4Evaporation;
  theMultiFragmentation = new G4StatMF;
  theFermiModel = new G4FermiBreakUp;
  thePhotonEvaporation = new G4PhotonEvaporation;
  SetParameters();
}

G4ExcitationHandler::~G4ExcitationHandler()
{
  if(MyOwnEvaporationClass) { delete theEvaporation; }
  delete theMultiFragmentation;
  delete theFermiModel;
  if(MyOwnPhotonEvaporationClass) { delete thePhotonEvaporation; }
}

void G4ExcitationHandler::SetParameters()
{
  //for inverse cross section choice
  theEvaporation->SetOPTxs(OPTxs);
  //for the choice of superimposed Coulomb Barrier for inverse cross sections
  theEvaporation->UseSICB(useSICB);
  theEvaporation->Initialise();
}

G4ReactionProductVector * G4ExcitationHandler::BreakItUp(const G4Fragment & theInitialState) const
{	
  
  // Variables existing until end of method
  G4Fragment * theInitialStatePtr = new G4Fragment(theInitialState);

  G4FragmentVector * theTempResult = 0;      // pointer which receives temporal results
  std::list<G4Fragment*> theEvapList;        // list to apply Evaporation or Fermi Break-Up
  std::list<G4bool> theFBUflag;              // list to apply Evaporation or Fermi Break-Up
  std::list<G4Fragment*> theEvapStableList;  // list to apply PhotonEvaporation
  std::list<G4Fragment*> theResults;         // list to store final result
  //
  //G4cout << "@@@@@@@@@@ Start G4Excitation Handler @@@@@@@@@@@@@" << G4endl;  
  //G4cout << theInitialState << G4endl;  
  
  // Variables to describe the excited configuration
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4int A = theInitialState.GetA_asInt();
  G4int Z = theInitialState.GetZ_asInt();

  G4NistManager* nist = G4NistManager::Instance();
  
  // JMQ 150909:  first step in de-excitation chain (SMM will be used only here)
  
  // In case A <= 1 the fragment will not perform any nucleon emission
  if (A <= 1)
    {
      theResults.push_back( theInitialStatePtr );
    }
  // check if a fragment is stable
  else if(exEnergy < minExcitation && nist->GetIsotopeAbundance(Z, A) > 0.0)
    {
      theResults.push_back( theInitialStatePtr );
    }  
  else  // In all cases apply once theFermiModel, theMultiFragmentation or theEvaporation
    {      
      // JMQ 150909: first step in de-excitation is treated separately 
      // Fragments after the first step are stored in theEvapList 
      // Statistical Multifragmentation will take place only once
      //
      // Fermi Break-Up always first
      G4bool wasFBU = false;
      if(A<maxAForFermiBreakUp && Z<maxZForFermiBreakUp) 
        {
          theTempResult = theFermiModel->BreakItUp(theInitialState);
	  // fragment was not decaied try to evaporate
          size_t nsec = theTempResult->size();
	  if(1 == nsec) {
            if((*theTempResult)[0] !=  theInitialStatePtr) { 
	      delete (*theTempResult)[0]; 
	    }
            delete theTempResult;
	    theTempResult = theEvaporation->BreakItUp(theInitialState);
	  } else if(nsec > 1) { wasFBU = true; }
        }
      else   if (exEnergy > minEForMultiFrag*A) 
        {
          theTempResult = theMultiFragmentation->BreakItUp(theInitialState);
	  // fragment was not decaied try to evaporate
	  if(1 == theTempResult->size()) {
            if((*theTempResult)[0] !=  theInitialStatePtr) { 
	      delete (*theTempResult)[0]; 
	    }
            delete theTempResult;
	    theTempResult = theEvaporation->BreakItUp(theInitialState);
	  }
        }
      else 
        {
          theTempResult = theEvaporation->BreakItUp(theInitialState);
        }
      
      G4bool deletePrimary = true;
      size_t nsec = theTempResult->size();
      if(nsec > 0) 
	{      
	  // Sort out secondary fragments
	  G4FragmentVector::iterator j;
	  for (j = theTempResult->begin(); j != theTempResult->end(); ++j)
	    {  
	      if((*j) == theInitialStatePtr) { deletePrimary = false; }
	      A = (*j)->GetA_asInt();  

	      if(A <= 1) { theResults.push_back(*j); }   // gamma, p, n
	      else if(1 == nsec) { theEvapStableList.push_back(*j); } // evaporation is not possible  
	      else  // Analyse fragment
		{
		  G4double exEnergy = (*j)->GetExcitationEnergy();
		  if(exEnergy < minExcitation) {
		    Z = (*j)->GetZ_asInt(); 
                    if(nist->GetIsotopeAbundance(Z, A) > 0.0) { 
		      theResults.push_back(*j); // stable fragment 
		    } else {   
		      theEvapList.push_back(*j); // unstable cold fragment
                      theFBUflag.push_back(wasFBU); 
		    }
		  } else { 
		    theEvapList.push_back(*j); // hot fragment 
		    theFBUflag.push_back(wasFBU); 
		  } 
		}
	    }
	}
      if( deletePrimary ) { delete theInitialStatePtr; }
      delete theTempResult;
    }
  //
  // JMQ 150909: Further steps in de-excitation chain follow ..
      
  //G4cout << "## After first step " << theEvapList.size() << " for evap;  "
  //	 << theEvapStableList.size() << " for photo-evap; " 
  //	 << theResults.size() << " results. " << G4endl; 

  // ------------------------------
  // De-excitation loop
  // ------------------------------
      
  std::list<G4Fragment*>::iterator iList;
  std::list<G4bool>::iterator ib = theFBUflag.begin();
  for (iList = theEvapList.begin(); iList != theEvapList.end(); ++iList, ++ib)
    {
      A = (*iList)->GetA_asInt(); 
      Z = (*iList)->GetZ_asInt();
	  
      // Fermi Break-Up 
      G4bool wasFBU = false;
      if ( !(*ib) && A < maxAForFermiBreakUp && Z < maxZForFermiBreakUp) 
	{
	  theTempResult = theFermiModel->BreakItUp(*(*iList));
	  // fragment was not decaied try to evaporate
          size_t nsec = theTempResult->size(); 
	  if(1 == nsec) {
            if((*theTempResult)[0] !=  theInitialStatePtr) { delete (*theTempResult)[0]; }
            delete theTempResult;
	    theTempResult = theEvaporation->BreakItUp(*(*iList));
	  } else if(nsec > 1) { wasFBU = true; }
	}
      else // apply Evaporation in another case
	{
	  theTempResult = theEvaporation->BreakItUp(*(*iList));
	}
      
      G4bool deletePrimary = true;
      size_t nsec = theTempResult->size();
		  
      // Sort out secondary fragments
      if ( nsec > 0 )
	{
	  G4FragmentVector::iterator j;
	  for (j = theTempResult->begin(); j != theTempResult->end(); ++j)
	    {
	      if((*j) == (*iList)) { deletePrimary = false; }
	      A = (*j)->GetA_asInt();

	      if(A <= 1) { theResults.push_back(*j); }                // gamma, p, n
	      else if(1 == nsec) { theEvapStableList.push_back(*j); } // no evaporation 
	      else  // Analyse fragment
		{
		  G4double exEnergy = (*j)->GetExcitationEnergy();
		  if(exEnergy < minExcitation) {
		    Z = (*j)->GetZ_asInt();
                    if(nist->GetIsotopeAbundance(Z, A) > 0.0) { 
		      theResults.push_back(*j); // stable fragment 
		    } else {   
		      theEvapList.push_back(*j); // unstable cold fragment 
                      theFBUflag.push_back(wasFBU); 
		    }
		  } else { 
		    theEvapList.push_back(*j); // hot fragment 
		    theFBUflag.push_back(wasFBU); 
		  }
		}
	    }
	}
      if( deletePrimary ) { delete (*iList); }
      delete theTempResult;
    } // end of the loop over theEvapList

  //G4cout << "## After 2nd step " << theEvapList.size() << " was evap;  "
  //	 << theEvapStableList.size() << " for photo-evap; " 
  //	 << theResults.size() << " results. " << G4endl; 
      
  // -----------------------
  // Photon-Evaporation loop
  // -----------------------
  
  // normally should not reach this point - it is kind of work around	      
  for(iList = theEvapStableList.begin(); iList != theEvapStableList.end(); ++iList)
    {
      // photon-evaporation is applied
      theTempResult = thePhotonEvaporation->BreakUpFragment(*iList);	  
      size_t nsec = theTempResult->size();
	  
      // if there is a gamma emission then
      if (nsec > 0)
	{
	  G4FragmentVector::iterator j;
	  for (j = theTempResult->begin(); j != theTempResult->end(); ++j)
	    {
	      theResults.push_back(*j); 
	    }
	}
      delete theTempResult;

      // priamry fragment is kept
      theResults.push_back(*iList); 

    } // end of photon-evaporation loop

  //G4cout << "## After 3d step " << theEvapList.size() << " was evap;  "
  //	 << theEvapStableList.size() << " was photo-evap; " 
  //	 << theResults.size() << " results. " << G4endl; 
    
  G4ReactionProductVector * theReactionProductVector = new G4ReactionProductVector;

  // MAC (24/07/08)
  // To optimise the storing speed, we reserve space in memory for the vector
  theReactionProductVector->reserve( theResults.size() );

  G4int theFragmentA, theFragmentZ;

  std::list<G4Fragment*>::iterator i;
  for (i = theResults.begin(); i != theResults.end(); ++i) 
    {
      theFragmentA = (*i)->GetA_asInt();
      theFragmentZ = (*i)->GetZ_asInt();
      G4ParticleDefinition* theKindOfFragment = 0;
      if (theFragmentA == 0) {       // photon or e-
	theKindOfFragment = (*i)->GetParticleDefinition();   
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
	theKindOfFragment = 
	  theTableOfIons->GetIon(theFragmentZ,theFragmentA,0.0);
      }
      if (theKindOfFragment != 0) 
	{
	  G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
	  theNew->SetMomentum((*i)->GetMomentum().vect());
	  theNew->SetTotalEnergy((*i)->GetMomentum().e());
	  theNew->SetFormationTime((*i)->GetCreationTime());
	  theReactionProductVector->push_back(theNew);
	}
      delete (*i);
    }

  return theReactionProductVector;
}

void G4ExcitationHandler::SetEvaporation(G4VEvaporation *const  value)
{
  delete theEvaporation;
  theEvaporation = value;
  MyOwnEvaporationClass = false;
  SetParameters();
}

void 
G4ExcitationHandler::SetMultiFragmentation(G4VMultiFragmentation *const  value)
{
  delete theMultiFragmentation;
  theMultiFragmentation = value;
}

void G4ExcitationHandler::SetFermiModel(G4VFermiBreakUp *const  value)
{
  delete theFermiModel;
  theFermiModel = value;
}

void 
G4ExcitationHandler::SetPhotonEvaporation(G4VEvaporationChannel *const  value)
{
  delete thePhotonEvaporation;
  thePhotonEvaporation = value;
  MyOwnPhotonEvaporationClass = false;
}

void G4ExcitationHandler::SetMaxZForFermiBreakUp(G4int aZ)
{
  maxZForFermiBreakUp = aZ;
}

void G4ExcitationHandler::SetMaxAForFermiBreakUp(G4int anA)
{
  maxAForFermiBreakUp = std::min(5,anA);
}

void G4ExcitationHandler::SetMaxAandZForFermiBreakUp(G4int anA, G4int aZ)
{
  SetMaxAForFermiBreakUp(anA);
  SetMaxZForFermiBreakUp(aZ);
}

void G4ExcitationHandler::SetMinEForMultiFrag(G4double anE)
{
  minEForMultiFrag = anE;
}






