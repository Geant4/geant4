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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
//
// Modified:
//
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
//      -according to Igor Pshenichnov, SMM will be applied (just in case) only once .
// (23 November 2009) By V.Ivanchenko:
//      - general cleanup of work with intermediate fragments
//

#include "G4ExcitationHandler.hh"
#include <list>

//#define debugphoton

G4ExcitationHandler::G4ExcitationHandler():
  // JMQ 160909 Fermi BreakUp & MultiFrag are on by default 
  // This is needed for activation of such models when G4BinaryLightIonReaction is used
  //  since no interface (for external activation via macro input file) is still available in this case.
  //maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),minEForMultiFrag(3.0*MeV),
  maxZForFermiBreakUp(1),maxAForFermiBreakUp(1),minEForMultiFrag(4.0*GeV),
  MyOwnEvaporationClass(true), MyOwnMultiFragmentationClass(true),MyOwnFermiBreakUpClass(true),
  MyOwnPhotonEvaporationClass(true),OPTxs(3),useSICB(false)
{                                                                          
  theTableOfParticles = G4ParticleTable::GetParticleTable();
  
  theEvaporation = new G4Evaporation;
  theMultiFragmentation = new G4StatMF;
  theFermiModel = new G4FermiBreakUp;
  thePhotonEvaporation = new G4PhotonEvaporation;
}

G4ExcitationHandler::G4ExcitationHandler(const G4ExcitationHandler &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4ExcitationHandler::copy_constructor: is meant to not be accessable! ");
}


G4ExcitationHandler::~G4ExcitationHandler()
{
  if (MyOwnEvaporationClass) delete theEvaporation;
  if (MyOwnMultiFragmentationClass) delete theMultiFragmentation;
  if (MyOwnFermiBreakUpClass) delete theFermiModel;
  if (MyOwnPhotonEvaporationClass) delete thePhotonEvaporation;
}

const G4ExcitationHandler & G4ExcitationHandler::operator=(const G4ExcitationHandler &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4ExcitationHandler::operator=: is meant to not be accessable! ");
  
  return *this;
}

G4bool G4ExcitationHandler::operator==(const G4ExcitationHandler &) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4ExcitationHandler::operator==: is meant to not be accessable! ");
  return false;
} 

G4bool G4ExcitationHandler::operator!=(const G4ExcitationHandler &) const
{
  throw G4HadronicException(__FILE__, __LINE__, "G4ExcitationHandler::operator!=: is meant to not be accessable! ");
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////
/// 25/07/08 16:45  Proposed by MAC ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector * G4ExcitationHandler::BreakItUp(const G4Fragment & theInitialState) const
{	
  //for inverse cross section choice
  theEvaporation->SetOPTxs(OPTxs);
  //for the choice of superimposed Coulomb Barrier for inverse cross sections
  theEvaporation->UseSICB(useSICB);
  
  // Variables existing until end of method
  //G4Fragment * theInitialStatePtr = const_cast<G4Fragment*>(&theInitialState);
  G4Fragment * theInitialStatePtr = new G4Fragment(theInitialState);

  G4FragmentVector * theTempResult = 0;      // pointer which receives temporal results
  std::list<G4Fragment*> theEvapList;        // list to apply Evaporation, SMF or Fermi Break-Up
  std::list<G4Fragment*> theEvapStableList;  // list to apply PhotonEvaporation
  std::list<G4Fragment*> theResults;         // list of final results

  std::list<G4Fragment*>::iterator iList;
  
  // Variables to describe the excited configuration
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4int A = static_cast<G4int>( theInitialState.GetA() +0.5 );
  G4int Z = static_cast<G4int>( theInitialState.GetZ() +0.5 );
  
  // JMQ 150909:  first step in de-excitation chain (SMM will be used only here)  
  // In case A <= 4 the fragment will not perform any nucleon emission
  if (A <= 4)
    {
      // I store G4Fragment* in theEvapStableList to apply thePhotonEvaporation later      
      theEvapStableList.push_back( theInitialStatePtr );
    }
  else  // If A > 4 we try to apply theFermiModel, theMultiFragmentation or theEvaporation
    {
      
      // JMQ 150909: first step in de-excitation is treated separately 
      // Fragments after the first step are stored in theEvapList 
      // Statistical Multifragmentation will take place (just in case) only here
      //
      // Test applicability
      // Initial State De-Excitation 
      if( A < GetMaxA() && Z < GetMaxZ() ) 
        {
          theTempResult = theFermiModel->BreakItUp(theInitialState);
        }
      else   if (exEnergy>GetMinE()*A) 
        {
          theTempResult = theMultiFragmentation->BreakItUp(theInitialState);
        }
      else 
        {
          theTempResult = theEvaporation->BreakItUp(theInitialState);
        }
      
      // Classify fragments
      G4FragmentVector::iterator j;
      for (j = theTempResult->begin(); j != theTempResult->end(); ++j)
	{  
	  if( (*j)->GetExcitationEnergy() > CLHEP::eV 
	      && ((*j)->GetA() > 4.5) )
	    { 	  
	      theEvapList.push_back(*j);
	    }
	  else
	    {
	      theEvapStableList.push_back(*j);
	    }
	}
      delete theTempResult;

      //
      // JMQ 150909: Further steps in de-excitation chain follow ..
      
      // ------------------------------
      // De-excitation loop
      // ------------------------------
      
      for (iList = theEvapList.begin(); iList != theEvapList.end(); ++iList)
	{
	  A = static_cast<G4int>((*iList)->GetA()+0.5);  // +0.5 to avoid bad truncation
	  Z = static_cast<G4int>((*iList)->GetZ()+0.5);
	  		  
	  if ( A < GetMaxA() && Z < GetMaxZ() ) // if satisfied apply Fermi Break-Up
	    {
	      theTempResult = theFermiModel->BreakItUp(*(*iList));
	    }
	  else // apply Evaporation in another case
	    {
	      theTempResult = theEvaporation->BreakItUp(*(*iList));
	    }
		
	  // one fragment in final state means no evaporation
	  size_t nevap = theTempResult->size();

	  // Classify fragments
	  for (j = theTempResult->begin(); j != theTempResult->end(); ++j)
	    {  
	      if ( (*j)->GetExcitationEnergy() > CLHEP::eV 
		  && (*j)->GetA() > 4.5 && 1 < nevap)
		{ 	  
		  theEvapList.push_back(*j);
		}
	      else
		{
		  theEvapStableList.push_back(*j);
		}
	    }

	  delete theTempResult;
	  delete (*iList); 
  
	} // end of the loop over theEvapList
    }
  
  // Now we try to deexcite by means of PhotonEvaporation those fragments
  // which are still excited.
  
  
  // -----------------------
  // Photon-Evaporation loop
  // -----------------------
  
  for (iList = theEvapStableList.begin(); iList != theEvapStableList.end(); ++iList)
    {
      A = static_cast<G4int>((*iList)->GetA()+0.5);
      exEnergy = (*iList)->GetExcitationEnergy();

      // stable fragment
      if ( (*iList)->GetA() < 1.5 || 
	   (*iList)->GetExcitationEnergy() < 1.1*CLHEP::eV ) 
	{
	  theResults.push_back(*iList);
	}
      else
	{
	  // photon-evaporation is applied
	  theTempResult = thePhotonEvaporation->BreakItUp(*(*iList));

	  // one fragment in final state means no evaporation
	  size_t nevap = theTempResult->size();

	  // Classify fragments
	  G4FragmentVector::iterator j;
	  for (j = theTempResult->begin(); j != theTempResult->end(); ++j)
	    {  
	      if ((*j)->GetExcitationEnergy() > CLHEP::eV 
		  && (*j)->GetA() > 4.5 && 1 < nevap)
		{ 	  
		  theEvapStableList.push_back(*j);
		}
	      else
		{
		  theResults.push_back(*j);
		}
	    }

	  delete theTempResult;
	  delete (*iList); 
	} // end of photon-evaporation loop
    } // end of loop on theEvapStableList
        
  // -----------------------
  // Fill final result
  // -----------------------

  // Change G4FragmentVector* to G4ReactionProductVector*
  G4ReactionProductVector * theReactionProductVector = 
    new G4ReactionProductVector();

  theReactionProductVector->reserve( theResults.size() );
    
  G4int idx = 0;
  for (iList = theResults.begin(); iList != theResults.end(); ++iList)
    {
      G4cout << "idx= " << idx << "  " << (*iList) << G4endl;
      G4int theFragmentA = static_cast<G4int>((*iList)->GetA());
      G4int theFragmentZ = static_cast<G4int>((*iList)->GetZ());
      G4LorentzVector theFragmentMomentum = (*iList)->GetMomentum();
      G4ParticleDefinition * theKindOfFragment = 0;
      if (theFragmentA == 0 && theFragmentZ == 0) {       // photon
	theKindOfFragment = G4Gamma::GammaDefinition();      
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
	theKindOfFragment = G4Alpha::AlphaDefinition();
      } else {
	theKindOfFragment = 
	  theTableOfParticles->FindIon(theFragmentZ,theFragmentA,0,theFragmentZ);
      }
      if (theKindOfFragment != 0) 
	{
	  G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
	  theNew->SetMomentum(theFragmentMomentum.vect());
	  theNew->SetTotalEnergy(theFragmentMomentum.e());
	  theNew->SetFormationTime((*iList)->GetCreationTime());
#ifdef PRECOMPOUND_TEST
	  theNew->SetCreatorModel((*iList)->GetCreatorModel());
#endif
	  theReactionProductVector->push_back(theNew);
	}
      delete (*iList); 
      ++idx; 
    }

#ifdef debug
  CheckConservation(theInitialState,theReactionProductVector);
#endif

  return theReactionProductVector;
}


#ifdef debug
void G4ExcitationHandler::CheckConservation(const G4Fragment & theInitialState,
					    G4FragmentVector * Result) const
{
  G4double ProductsEnergy =0;
  G4ThreeVector ProductsMomentum;
  G4int ProductsA = 0;
  G4int ProductsZ = 0;
  G4FragmentVector::iterator h;
  for (h = Result->begin(); h != Result->end(); h++) {
    G4LorentzVector tmp = (*h)->GetMomentum();
    ProductsEnergy += tmp.e();
    ProductsMomentum += tmp.vect();
    ProductsA += static_cast<G4int>((*h)->GetA());
    ProductsZ += static_cast<G4int>((*h)->GetZ());
  }
  
  if (ProductsA != theInitialState.GetA()) {
    G4cout << "!!!!!!!!!! Baryonic Number Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4ExcitationHandler.cc: Barionic Number Conservation test for deexcitation fragments" 
	   << G4endl; 
    G4cout << "Initial A = " << theInitialState.GetA() 
	   << "   Fragments A = " << ProductsA << "   Diference --> " 
	   << theInitialState.GetA() - ProductsA << G4endl;
  }
  if (ProductsZ != theInitialState.GetZ()) {
    G4cout << "!!!!!!!!!! Charge Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4ExcitationHandler.cc: Charge Conservation test for deexcitation fragments" 
	   << G4endl; 
    G4cout << "Initial Z = " << theInitialState.GetZ() 
	   << "   Fragments Z = " << ProductsZ << "   Diference --> " 
	   << theInitialState.GetZ() - ProductsZ << G4endl;
  }
  if (std::abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4ExcitationHandler.cc: Energy Conservation test for deexcitation fragments" 
	   << G4endl; 
    G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	   << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	   << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
  } 
  if (std::abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
      std::abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
      std::abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Momentum Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4ExcitationHandler.cc: Momentum Conservation test for deexcitation fragments" 
	   << G4endl; 
    G4cout << "Initial P = " << theInitialState.GetMomentum().vect() << " MeV"
	   << "   Fragments P = " << ProductsMomentum  << " MeV   Diference --> " 
	   << theInitialState.GetMomentum().vect() - ProductsMomentum << " MeV" << G4endl;
  }
  return;
}
#endif




