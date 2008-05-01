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
// Modif (25 Abr 2007) by M. A. Cortes Giraldo:
//      -Max Z,A for Fermi Break-Up turns to 9,17 by default
//      -Bug fixed on BreakItUp()
// Modif (30 June 1998) by V. Lara:
//      -Modified the Transform method for use G4ParticleTable and 
//       therefore G4IonTable. It makes possible to convert all kind 
//       of fragments (G4Fragment) produced in deexcitation to 
//       G4DynamicParticle
//      -It uses default algorithms for:
//              Evaporation: G4Evaporation
//              MultiFragmentation: G4StatMF 
//              Fermi Breakup model: G4FermiBreakUp


#include "G4ExcitationHandler.hh"
#include <list>

//#define debugphoton


G4ExcitationHandler::G4ExcitationHandler():
  // Fermi BreakUp is on and MultiFrag is off by default
  maxZForFermiBreakUp(9),maxAForFermiBreakUp(17),minEForMultiFrag(4.0*GeV),
  MyOwnEvaporationClass(true), MyOwnMultiFragmentationClass(true),MyOwnFermiBreakUpClass(true),
  MyOwnPhotonEvaporationClass(true)
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

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// 28/01/08 13:30  Proposed by MAC /////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

G4ReactionProductVector * G4ExcitationHandler::BreakItUp(const G4Fragment & theInitialState) const
{

  G4FragmentVector * theResult = 0;         // to be returned (will be created dynamically later)

  // Variables to describe the initial state
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4int A = static_cast<G4int>(theInitialState.GetA()+0.5);  // +0.5 to avoid bad truncation
  G4int Z = static_cast<G4int>(theInitialState.GetZ()+0.5);

  // Variables existing until end of method
  G4Fragment theExcitedNucleus;  // object to be passed in BreakItUp methods

  // In case A <= 4 the fragment will not perform particle emission
  if (A <= 4)
    {
      theResult = new G4FragmentVector();
      theResult->push_back( new G4Fragment(theInitialState) );
    }
  else  // If A > 4 we can apply this algorithm
    {
      // Initial State De-excitation
      if ( A < GetMaxA() && Z < GetMaxZ() )  // if satisfied we apply Fermi Break-Up
	{
	  theResult = theFermiModel->BreakItUp(theInitialState);
	}
      else if (exEnergy > GetMinE()*A)   // if satisfied we apply MF instead of Evaporation
	{
	  theResult = theMultiFragmentation->BreakItUp(theInitialState);
	}
      else   // Evaporation model (there is competition between Evaporation, Fission, PhotonEvaporation!)
	{
	  theResult = theEvaporation->BreakItUp(theInitialState);
	}

      // De-excitation loop
      // ------------------

      G4FragmentVector * theTempResult = 0; // pointer which receives temporal results
      for (G4FragmentVector::iterator i = theResult->begin(); i != theResult->end(); i++)
	{
	  exEnergy = (*i)->GetExcitationEnergy();
	  if (exEnergy > 0.0)
	    {
	      // These commented lines are to be used if we make conditions about them
	      // A = static_cast<G4int>((*i)->GetA()+0.5);
	      // Z = static_cast<G4int>((*i)->GetZ()+0.5);
	      theExcitedNucleus = *(*i); // G4Fragment type needed to pass object in BreakItUp method

	      // only valid if using evaporation model for the fragments is enough
	      theTempResult = theEvaporation->BreakItUp(theExcitedNucleus);

	      if (theTempResult->size() > 1)
		// if so, some deexcitation has been made
		{
		  // first we add final state to theResult vector
		  for (G4FragmentVector::reverse_iterator ri = theTempResult->rbegin(); 
		       ri != theTempResult->rend(); ++ri)
		    // reverse is applied in order to have residual nucleus in first position
		    {
		      theResult->push_back(*ri);
		    }
		  // now we remove initial state
		  delete (*i);     // in the dynamic memory
		  i = theResult->erase(i); // and now erase the pointer in the vector
		  i--;    // if we don't do this, the following object will be skipped
		  
		  delete theTempResult;
		}
	      else   // if not, simply I have to clean up the temporal vector
		{
		  for (G4FragmentVector::iterator j = theTempResult->begin(); 
		       j != theTempResult->end(); j++)
		    {delete (*j);}
		  delete theTempResult;
		}
	    } // end if (exEnergy > 0) block
	  // increment iterator if exEnergy <= 0

	} // END OF LOOP

    } // if/else A<=4 block

  // Now we try to deexcite by means of PhotonEvaporation those fragments
  // which are excited.

  G4FragmentVector* theTempResult = 0;

  // Photon-Evaporation loop
  // ----------------------
  for (G4FragmentVector::iterator j = theResult->begin(); j != theResult->end(); j++)
    {
      A = static_cast<G4int>((*j)->GetA()+0.5);
      exEnergy = (*j)->GetExcitationEnergy();
      if ( A > 1 && exEnergy > 0.1*eV ) // if so, photon-evaporation is applied
	{
	  theExcitedNucleus = *(*j);
	  theTempResult = thePhotonEvaporation->BreakItUp(theExcitedNucleus);

	  // if Gamma Evaporation has succeed then
	  if (theTempResult->size() > 1)
	    {
	      // first we add final state to theResult vector
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

		  theResult->push_back(*ri);
		}
	     
	      // now we remove initial state
	      delete (*j);     // in dynamic memory
	      j = theResult->erase(j); // and now erase the pointer in the vector
	      j--;    // if we don't do this, the following object will be skipped
	      delete theTempResult;
	    }
	  else         // if not succeed, we deexcited manually with one gamma
	    {
	      for (G4FragmentVector::iterator i = theTempResult->begin(); i != theTempResult->end(); i++)
		{
		  delete (*i);
		}
	      delete theTempResult;

#ifdef debugphoton
              G4cout << "G4ExcitationHandler: Gamma Evaporation could not deexcite the nucleus: \n"
                     << "-----------------------------------------------------------------------\n"
                     << theExcitedNucleus << '\n'
                     << "-----------------------------------------------------------------------\n";
#endif

	      // Let's create a G4Fragment pointer representing the gamma emmited
	      G4double GammaEnergy = (*j)->GetExcitationEnergy();
	      G4double cosTheta = 1. - 2. * G4UniformRand();
	      G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
              G4double phi = twopi * G4UniformRand();
              G4ThreeVector GammaP(GammaEnergy * sinTheta * std::cos(phi),
                                   GammaEnergy * sinTheta * std::sin(phi),
                                   GammaEnergy * cosTheta );
              G4LorentzVector Gamma4P(GammaP,GammaEnergy);
              G4Fragment * theHandlerPhoton = new G4Fragment(Gamma4P,G4Gamma::GammaDefinition());

	      // And now we update momentum and energy for the nucleus
              G4double Mass = (*j)->GetGroundStateMass();
              G4ThreeVector ResidualP((*j)->GetMomentum().vect() - GammaP);
              G4double ResidualE = std::sqrt(ResidualP*ResidualP + Mass*Mass);
              G4LorentzVector Residual4P(ResidualP,ResidualE);
              (*j)->SetMomentum(Residual4P); // Now this fragment has been deexcited

#ifdef PRECOMPOUND_TEST
	      theHandlerPhoton->SetCreatorModel("G4ExcitationHandler");
#endif

	      // Finally, we add the gamma to the list
	      theResult->push_back(theHandlerPhoton); 

#ifdef debugphoton
              G4cout << "Emmited photon:\n"
                     << theResult->back() << '\n'
                     << "Residual nucleus after photon emission:\n"
                     << *(*j) << '\n'
                     << "-----------------------------------------------------------------------\n";
#endif

	    }
	}
      // increment iterator if a proton, gamma or excitation energy small
      // so we don't have to do anything 

    } // end of photon-evaporation loop

#ifdef debug
  CheckConservation(theInitialState,theResult);
#endif

  // Change G4FragmentVector* by G4ReactionProductVector*
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




