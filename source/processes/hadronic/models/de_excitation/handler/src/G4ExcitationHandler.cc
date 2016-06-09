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
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option
//
// Modif (24 Jul 2008) by M. A. Cortes Giraldo:
//      -Max Z,A for Fermi Break-Up turns to 9,17 by default
//      -BreakItUp() reorganised and bug in Evaporation loop fixed
//      -Transform() optimised
// Modif (30 June 1998) by V. Lara:
//      -Modified the Transform method for use G4ParticleTable and 
//       therefore G4IonTable. It makes possible to convert all kind 
//       of fragments (G4Fragment) produced in deexcitation to 
//       G4DynamicParticle
//      -It uses default algorithms for:
//              Evaporation: G4Evaporation
//              MultiFragmentation: G4StatMF 
//              Fermi Breakup model: G4FermiBreakUp
//
// Modif (03 September 2008) by J. M. Quesada for external choice of inverse 
// cross section option (default OPTxs=3)
// JMQ (06 September 2008) Also external choices have been added for 
// superimposed Coulomb barrier (if useSICBis set true, by default is false) 

#include "G4ExcitationHandler.hh"
#include <list>

//#define debugphoton


G4ExcitationHandler::G4ExcitationHandler():
  // Fermi BreakUp is on and MultiFrag is off by default
  //  maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),minEForMultiFrag(4.0*GeV),
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

  // Pointer which will be used to return the final production vector
  G4FragmentVector * theResult = 0;

  // Variables existing until end of method
  G4Fragment * theInitialStatePtr = new G4Fragment(theInitialState);
  G4Fragment theExcitedNucleus;              // object to be passed in BreakItUp methods
  G4FragmentVector * theTempResult = 0;      // pointer which receives temporal results
  std::list<G4Fragment*> theEvapList;        // list to apply Evaporation, SMF or Fermi Break-Up
  std::list<G4Fragment*> theEvapStableList;  // list to apply PhotonEvaporation
  std::list<G4Fragment*> theFinalStableList; // list to store final result
  std::list<G4Fragment*>::iterator iList;

  // Variables to describe the excited configuration
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4int A = static_cast<G4int>( theInitialState.GetA() +0.5 );
  G4int Z = static_cast<G4int>( theInitialState.GetZ() +0.5 );

  // In case A <= 4 the fragment will not perform any nucleon emission
  if (A <= 4)
    {
      // I store G4Fragment* in theEvapStableList to apply thePhotonEvaporation later

      theEvapStableList.push_back( theInitialStatePtr );
    }

  else  // If A > 4 we try to apply theFermiModel, theMultiFragmentation or theEvaporation
    {

      // Store original state in theEvapList
      theEvapList.push_back( theInitialStatePtr );

      // ------------------------------
      // De-excitation loop
      // ------------------------------

      for (iList = theEvapList.begin(); iList != theEvapList.end(); iList++)
	{
	  exEnergy = (*iList)->GetExcitationEnergy();

	  if (exEnergy > 0.0)
	    {
	      // Check conditions for each model
	      theExcitedNucleus = *(*iList);
	      A = static_cast<G4int>((*iList)->GetA()+0.5);  // +0.5 to avoid bad truncation
	      Z = static_cast<G4int>((*iList)->GetZ()+0.5);

	      if ( A < GetMaxA() && Z < GetMaxZ() ) // if satisfied apply Fermi Break-Up
		{
		  theTempResult = theFermiModel->BreakItUp(theExcitedNucleus);
		}
	      else if (exEnergy > GetMinE()*A)   // if satisfied apply SMF
		{
		  theTempResult = theMultiFragmentation->BreakItUp(theExcitedNucleus);
		}
	      else // apply Evaporation in another case
		{
		  theTempResult = theEvaporation->BreakItUp(theExcitedNucleus);
		}

	      // New configuration is stored in theTempResult, so we can free
	      // the memory where the previous configuration is

	      delete (*iList);

	      // And now the theTempResult->size() tells us if the configuration has changed

	      if ( theTempResult->size() > 1 )
		{
		  // push_back the result to the end of theEvapList (this same list)

		  for (G4FragmentVector::iterator j = theTempResult->begin();
		       j != theTempResult->end(); j++)
		    {
		      theEvapList.push_back(*j);
		    }
		}
	      else
		{
		  // push_back the result to theEvapStableList, because
		  // is still excited, but cannot emmit more nucleons

		  for (G4FragmentVector::iterator j = theTempResult->begin();
		       j != theTempResult->end(); j++)
		    {
		      theEvapStableList.push_back(*j);
		    }
		}

	      // after working with theTempResult, clear and delete it
	      theTempResult->clear();
	      delete theTempResult;

	    }
	  else // exEnergy = 0.0
	    {
	      // if this fragment is at ground state,
	      // store it in theFinalStableList

	      theFinalStableList.push_back(*iList);

	    } // endif (exEnergy > 0.0)

	} // end of the loop over theEvapList

      theEvapList.clear();   // clear all the list and do not free memory pointed by
                             // each element because this have been done before!
    }

  // Now we try to deexcite by means of PhotonEvaporation those fragments
  // which are still excited.


  // -----------------------
  // Photon-Evaporation loop
  // -----------------------

  for (iList = theEvapStableList.begin(); iList != theEvapStableList.end(); iList++)
    {
      A = static_cast<G4int>((*iList)->GetA()+0.5);
      exEnergy = (*iList)->GetExcitationEnergy();

      if ( A > 1 && exEnergy > 0.1*eV ) // if so, photon-evaporation is applied
	{
	  theExcitedNucleus = *(*iList);
	  theTempResult = thePhotonEvaporation->BreakItUp(theExcitedNucleus);

	  // if there is a gamma emission then
	  if (theTempResult->size() > 1)
	    {
	      // first free the memory occupied by the previous state
	      delete (*iList);

	      // and now add the final state from gamma emission to the end of
	      // theEvapStableList
	      for (G4FragmentVector::reverse_iterator ri = theTempResult->rbegin(); 
		   ri != theTempResult->rend(); ++ri)
		// reversed is applied in order to have residual nucleus in first position
		{

#ifdef PRECOMPOUND_TEST
		  if ((*ri)->GetA() == 0)
		    (*ri)->SetCreatorModel(G4String("G4PhotonEvaporation"));
		  else
		    (*ri)->SetCreatorModel(G4String("ResidualNucleus"));
#endif

		  theEvapStableList.push_back(*ri);
		}
	     
	      // now we clean and remove the temporal vector
	      theTempResult->clear();
	      delete theTempResult;

	    }
	  else  // if theTempResult->size() = 1
	    {
	      // if there is not any gamma emission from this excited fragment
	      // we have to emmit a gamma which forces the deexcitation

	      // First I clean completely theTempResult

	      for (G4FragmentVector::iterator j = theTempResult->begin();
		   j != theTempResult->end(); j++)
		{
		  delete (*j);
		}

	      theTempResult->clear();
	      delete theTempResult;

#ifdef debugphoton
              G4cout << "G4ExcitationHandler: Gamma Evaporation could not deexcite the nucleus: \n"
                     << "-----------------------------------------------------------------------\n"
                     << theExcitedNucleus << '\n'
                     << "-----------------------------------------------------------------------\n";
#endif

	      // Let's create a G4Fragment pointer representing the gamma emmited
	      G4double GammaEnergy = (*iList)->GetExcitationEnergy();
	      G4double cosTheta = 1. - 2. * G4UniformRand();
	      G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
              G4double phi = twopi * G4UniformRand();
              G4ThreeVector GammaP(GammaEnergy * sinTheta * std::cos(phi),
                                   GammaEnergy * sinTheta * std::sin(phi),
                                   GammaEnergy * cosTheta );
              G4LorentzVector Gamma4P(GammaP,GammaEnergy);
              G4Fragment * theHandlerPhoton = new G4Fragment(Gamma4P,G4Gamma::GammaDefinition());

	      // And now we update momentum and energy for the nucleus
              G4double Mass = (*iList)->GetGroundStateMass();
              G4ThreeVector ResidualP((*iList)->GetMomentum().vect() - GammaP);
              G4double ResidualE = std::sqrt(ResidualP*ResidualP + Mass*Mass);
              G4LorentzVector Residual4P(ResidualP,ResidualE);
              (*iList)->SetMomentum(Residual4P); // Now this fragment has been deexcited!

	      // we store the deexcited fragment in theFinalStableList
	      theFinalStableList.push_back(*iList);

#ifdef PRECOMPOUND_TEST
	      theHandlerPhoton->SetCreatorModel("G4ExcitationHandler");
#endif

	      // Finally, we add theHandlerPhoton to theFinalStableList
	      theFinalStableList.push_back(theHandlerPhoton); 

#ifdef debugphoton
              G4cout << "Emmited photon:\n"
                     << theFinalStableList.back() << '\n'
                     << "Residual nucleus after photon emission:\n"
                     << *(*iList) << '\n'
                     << "-----------------------------------------------------------------------\n";
#endif

	    }
	}
      else // case of a nucleon, gamma or very small excitation energy
	{
	  // we don't have to do anything, just store the fragment in theFinalStableList

	  theFinalStableList.push_back(*iList);
	
	} // A > 1 && exEnergy > 0.1*eV

    } // end of photon-evaporation loop


  // The deexcitation from fragments inside theEvapStableList has been finished, so...
  theEvapStableList.clear();

  // Now the final state is in theFinalStableList, and we have to send it to theResult vector

  theResult = new G4FragmentVector;
  theResult->reserve( theFinalStableList.size() * sizeof(G4Fragment*) );
  // We reserve enough memory to optimise the storing speed

  for (iList = theFinalStableList.begin(); iList != theFinalStableList.end(); iList++)
    {
      theResult->push_back(*iList);
    }

  // After storing the final state , we can clear theFinalStableList
  theFinalStableList.clear();

#ifdef debug
  CheckConservation(theInitialState,theResult);
#endif

  // Change G4FragmentVector* to G4ReactionProductVector*
  return Transform(theResult);
}



G4ReactionProductVector * 
G4ExcitationHandler::Transform(G4FragmentVector * theFragmentVector) const
{
  if (theFragmentVector == 0) return 0;
  
  // Conversion from G4FragmentVector to G4ReactionProductVector
  G4ParticleDefinition *theGamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
  G4ParticleDefinition *theProton = G4Proton::ProtonDefinition();   
  G4ParticleDefinition *theDeuteron = G4Deuteron::DeuteronDefinition();
  G4ParticleDefinition *theTriton = G4Triton::TritonDefinition();
  G4ParticleDefinition *theHelium3 = G4He3::He3Definition();
  G4ParticleDefinition *theAlpha = G4Alpha::AlphaDefinition();
  G4ParticleDefinition *theKindOfFragment = 0;
  theNeutron->SetVerboseLevel(2);
  G4ReactionProductVector * theReactionProductVector = new G4ReactionProductVector;

  // MAC (24/07/08)
  // To optimise the storing speed, we reserve space in memory for the vector
  theReactionProductVector->reserve( theFragmentVector->size() * sizeof(G4ReactionProduct*) );

  G4int theFragmentA, theFragmentZ;
  G4LorentzVector theFragmentMomentum;

  G4FragmentVector::iterator i;
  for (i = theFragmentVector->begin(); i != theFragmentVector->end(); i++) {
    //    std::cout << (*i) <<'\n';
    theFragmentA = static_cast<G4int>((*i)->GetA());
    theFragmentZ = static_cast<G4int>((*i)->GetZ());
    theFragmentMomentum = (*i)->GetMomentum();
    theKindOfFragment = 0;
    if (theFragmentA == 0 && theFragmentZ == 0) {       // photon
      theKindOfFragment = theGamma;      
    } else if (theFragmentA == 1 && theFragmentZ == 0) { // neutron
      theKindOfFragment = theNeutron;
    } else if (theFragmentA == 1 && theFragmentZ == 1) { // proton
      theKindOfFragment = theProton;
    } else if (theFragmentA == 2 && theFragmentZ == 1) { // deuteron
      theKindOfFragment = theDeuteron;
    } else if (theFragmentA == 3 && theFragmentZ == 1) { // triton
      theKindOfFragment = theTriton;
    } else if (theFragmentA == 3 && theFragmentZ == 2) { // helium3
      theKindOfFragment = theHelium3;
    } else if (theFragmentA == 4 && theFragmentZ == 2) { // alpha
      theKindOfFragment = theAlpha;
    } else {
      theKindOfFragment = theTableOfParticles->FindIon(theFragmentZ,theFragmentA,0,theFragmentZ);
    }
    if (theKindOfFragment != 0) 
      {
	G4ReactionProduct * theNew = new G4ReactionProduct(theKindOfFragment);
	theNew->SetMomentum(theFragmentMomentum.vect());
	theNew->SetTotalEnergy(theFragmentMomentum.e());
	theNew->SetFormationTime((*i)->GetCreationTime());
#ifdef PRECOMPOUND_TEST
	theNew->SetCreatorModel((*i)->GetCreatorModel());
#endif
	theReactionProductVector->push_back(theNew);
      }
  }
  if (theFragmentVector != 0)
    { 
      std::for_each(theFragmentVector->begin(), theFragmentVector->end(), DeleteFragment());
      delete theFragmentVector;
    }
  G4ReactionProductVector::iterator debugit;
  for(debugit=theReactionProductVector->begin(); 
      debugit!=theReactionProductVector->end(); debugit++)
    {
    if((*debugit)->GetTotalEnergy()<1.*eV)
      {
	if(getenv("G4DebugPhotonevaporationData"))
	  {
	    G4cerr << "G4ExcitationHandler: Warning: Photonevaporation data not exact."<<G4endl;
	    G4cerr << "G4ExcitationHandler: Warning: Found gamma with energy = "
		   << (*debugit)->GetTotalEnergy()/MeV << "MeV"
		   << G4endl;
	  }
	delete (*debugit);
	*debugit = 0;
      }
  }
  G4ReactionProduct* tmpPtr=0;
  theReactionProductVector->erase(std::remove_if(theReactionProductVector->begin(),
                                                 theReactionProductVector->end(),
                                                 std::bind2nd(std::equal_to<G4ReactionProduct*>(),
                                                              tmpPtr)),
				  theReactionProductVector->end());
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




