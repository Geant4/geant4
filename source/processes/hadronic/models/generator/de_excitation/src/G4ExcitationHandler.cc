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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (May 1998)
// Modif (30 June 1998) by V. Lara:
//      -Modified the Transform method for use G4ParticleTable and 
//       therefore G4IonTable. It makes possible to convert all kind 
//       of fragments (G4Fragment) produced in deexcitation to 
//       G4DynamicParticle
//      -It uses default algorithms for:
//              Evaporation: G4StatEvaporation
//              MultiFragmentation: G4DummyMF (a dummy one)
//              Fermi Breakup model: G4StatFermiBreakUp


#include "G4ExcitationHandler.hh"
#include "g4std/list"

G4ExcitationHandler::G4ExcitationHandler():
  maxZForFermiBreakUp(8),maxAForFermiBreakUp(16),minEForMultiFrag(1000.0*MeV), // make Multifrag. unavailable 
  MyOwnEvaporationClass(true), MyOwnMultiFragmentationClass(true),MyOwnFermiBreakUpClass(true),
  MyOwnPhotonEvaporationClass(true)
{                                                                          // change by 3.0
  theTableOfParticles = G4ParticleTable::GetParticleTable();

  theEvaporation = new G4Evaporation;
  theMultiFragmentation = new G4StatMF;
  theFermiModel = new G4FermiBreakUp;
  thePhotonEvaporation = new G4PhotonEvaporation;
}

G4ExcitationHandler::G4ExcitationHandler(const G4ExcitationHandler &right)
{
  G4Exception("G4ExcitationHandler::copy_constructor: is meant to not be accessable! ");
}


G4ExcitationHandler::~G4ExcitationHandler()
{
  if (MyOwnEvaporationClass) delete theEvaporation;
  if (MyOwnMultiFragmentationClass) delete theMultiFragmentation;
  if (MyOwnFermiBreakUpClass) delete theFermiModel;
  if (MyOwnPhotonEvaporationClass) delete thePhotonEvaporation;
}


const G4ExcitationHandler & G4ExcitationHandler::operator=(const G4ExcitationHandler &right)
{
  G4Exception("G4ExcitationHandler::operator=: is meant to not be accessable! ");

  return *this;
}


G4bool G4ExcitationHandler::operator==(const G4ExcitationHandler &right) const
{
  G4Exception("G4ExcitationHandler::operator==: is meant to not be accessable! ");
  return false;
} 

G4bool G4ExcitationHandler::operator!=(const G4ExcitationHandler &right) const
{
  G4Exception("G4ExcitationHandler::operator!=: is meant to not be accessable! ");
  return true;
}


G4ReactionProductVector * G4ExcitationHandler::BreakItUp(const G4Fragment &theInitialState) const
{

  G4FragmentVector * theResult = 0; 
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4double A = theInitialState.GetA();
  G4int Z = G4int(theInitialState.GetZ());
  G4FragmentVector* theTempResult = 0; 
  G4Fragment theExcitedNucleus;

  // Test applicability
  if (A > 4) {
    // Initial State De-Excitation 
    if(A<GetMaxA()&&Z<GetMaxZ()&&
       exEnergy>G4NucleiPropertiesTable::GetBindingEnergy(Z,A)) {
    
      theResult = theFermiModel->BreakItUp(theInitialState);
    
    } else   if (exEnergy>GetMinE()*A) {
      
      theResult = theMultiFragmentation->BreakItUp(theInitialState);
    
    } else {
    
      theResult = theEvaporation->BreakItUp(theInitialState);
    }
    



    // De-Excitation loop
    // ------------------
    // Check if there are excited fragments
    G4std::list<G4Fragment*> theResultList;
    G4FragmentVector::iterator j;
    G4std::list<G4Fragment*>::iterator i;
    for (j = theResult->begin(); j != theResult->end();j++) 
      theResultList.push_back(*j);
    theResult->clear();
    for (i = theResultList.begin(); i != theResultList.end(); i++) {
      exEnergy = (*i)->GetExcitationEnergy();
      if (exEnergy > 0.0) {
	A = (*i)->GetA();
	Z = G4int((*i)->GetZ());
	theExcitedNucleus = *(*i);
	// try to de-excite this fragment
	if(A<GetMaxA()&&Z<GetMaxZ()&&
	   exEnergy>G4NucleiPropertiesTable::GetBindingEnergy(Z,A)) {
	  // Fermi Breakup
	  theTempResult = theFermiModel->BreakItUp(theExcitedNucleus);
	} else   if(exEnergy>GetMinE()*A) {
	  // Multifragmentation
	  theTempResult = theMultiFragmentation->BreakItUp(theExcitedNucleus);
	} else {
	  // Evaporation
	  theTempResult = theEvaporation->BreakItUp(theExcitedNucleus);
	}
	// The Nucleus has been fragmented?
	if (theTempResult->size() > 1) {
	  // If so :

	  // Remove excited fragment from the result 
	  //	delete theResult->removeAt(i--);
	  delete (*i);
	  i = theResultList.erase(i--);
	  // and add theTempResult elements to theResult
	  while (!theTempResult->empty()) {
	    theResultList.push_back(*(theTempResult->end()-1));
	    theTempResult->pop_back();
	  }
	  delete theTempResult;
	} else { // If not :
	  // it doesn't matter, we Follow with the next fragment but
	  // I have to make 
	  while ( !theTempResult->empty() ) {
	    delete *(theTempResult->end()-1);
	    theTempResult->pop_back();
	  }
	  delete theTempResult;
	}
      }
    }
    for (i = theResultList.begin(); i != theResultList.end(); i++)
      theResult->push_back(*i);
    theResultList.clear();
  } else {  // if A > 4
    theResult->push_back(new G4Fragment(theInitialState));
  }

  // Now we try to deexcite by means of PhotonEvaporation those fragments
  // which are excited.

  theTempResult = 0;
  G4std::list<G4Fragment*> theResultList;
  G4std::list<G4Fragment*>::iterator j;
  G4FragmentVector::iterator i;
  for (i = theResult->begin(); i != theResult->end();i++) 
    theResultList.push_back(*i);
  theResult->clear();
  
  for (j = theResultList.begin(); j != theResultList.end(); j++) {
    if ((*j)->GetExcitationEnergy() > 0.0 && (*j)->GetA() > 1) {
      theExcitedNucleus = *(*j);
      theTempResult = thePhotonEvaporation->BreakItUp(theExcitedNucleus);
      // If Gamma Evaporation has succeed then
      if (theTempResult->size() > 1) {
	// Remove excited fragment from the result 
	delete (*j);
        theResultList.erase(j--);
	// and add theTempResult elements to theResult
	while (!theTempResult->empty()) {
	  theResultList.push_back(*(theTempResult->end()-1));
	  theTempResult->pop_back();
	}
	delete theTempResult;
      }
      // In other case, just clean theTempResult and continue
      else {
	while (!theTempResult->empty()) {
	  delete *(theTempResult->end()-1);
	  theTempResult->pop_back();
	}
	delete theTempResult;
#ifdef debug
	G4cout << "G4ExcitationHandler: Gamma Evaporation could not deexcite the nucleus: "
	       << G4endl
	       << "-----------------------------------------------------------------------"  
	       << G4endl;
	G4cout << theExcitedNucleus
	       << G4endl
	       << "-----------------------------------------------------------------------"  
	       << G4endl;
#endif
	G4double GammaEnergy = (*j)->GetExcitationEnergy();
	G4double cosTheta = 1. - 2. * G4UniformRand();
	G4double sinTheta = sqrt(1. - cosTheta * cosTheta);
	G4double phi = twopi * G4UniformRand();
	G4ThreeVector GammaP(GammaEnergy * sinTheta * cos(phi),
			     GammaEnergy * sinTheta * sin(phi),
			     GammaEnergy * cosTheta );
	
	
	G4LorentzVector Gamma4P(GammaP,GammaEnergy);
	G4double Mass = (*j)->GetGroundStateMass();
	G4ThreeVector ResidualP = (*j)->GetMomentum().vect();
	G4LorentzVector Residual4P(ResidualP,sqrt(Mass*Mass+ResidualP.mag2()));
	(*j)->SetMomentum(Residual4P);
	theResultList.push_back( new G4Fragment(Gamma4P,G4Gamma::GammaDefinition()) );
      }	
    } 
  }
  for (j = theResultList.begin(); j != theResultList.end(); j++)
    theResult->push_back(*j);
  theResultList.clear();
  
  
#ifdef debug
  CheckConservation(theInitialState,theResult);
#endif
  // Change G4FragmentVector by G4DynamicParticle
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
      theFragmentA = G4int((*i)->GetA());
      theFragmentZ = G4int((*i)->GetZ());
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
#ifdef pctest
	  theNew->SetCreatorModel((*i)->GetCreatorModel());
#endif
	  theReactionProductVector->push_back(theNew);
      }
  }
  if (theFragmentVector != 0)
  { 
      while (!theFragmentVector->empty()) {
	  delete *(theFragmentVector->end()-1);
	  theFragmentVector->pop_back();
      }
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
	  theReactionProductVector->erase(debugit);
	  if(debugit!=theReactionProductVector->begin()) debugit--;
      }
  }
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
    ProductsA += G4int((*h)->GetA());
    ProductsZ += G4int((*h)->GetZ());
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
  if (abs(ProductsEnergy-theInitialState.GetMomentum().e()) > 1.0*keV) {
    G4cout << "!!!!!!!!!! Energy Conservation Violation !!!!!!!!!!" << G4endl;
    G4cout << "G4ExcitationHandler.cc: Energy Conservation test for deexcitation fragments" 
	   << G4endl; 
    G4cout << "Initial E = " << theInitialState.GetMomentum().e()/MeV << " MeV"
	   << "   Fragments E = " << ProductsEnergy/MeV  << " MeV   Diference --> " 
	   << (theInitialState.GetMomentum().e() - ProductsEnergy)/MeV << " MeV" << G4endl;
  } 
  if (abs(ProductsMomentum.x()-theInitialState.GetMomentum().x()) > 1.0*keV || 
      abs(ProductsMomentum.y()-theInitialState.GetMomentum().y()) > 1.0*keV ||
      abs(ProductsMomentum.z()-theInitialState.GetMomentum().z()) > 1.0*keV) {
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




