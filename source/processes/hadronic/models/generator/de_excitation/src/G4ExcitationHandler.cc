// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ExcitationHandler.cc,v 1.3 1999/05/28 17:14:41 hpw Exp $
// GEANT4 tag $Name: geant4-00-01 $
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

G4ExcitationHandler::G4ExcitationHandler():MyOwnEvaporationClass(true),
  MyOwnMultiFragmentationClass(true),MyOwnFermiBreakUpClass(true),
  MyOwnPhotonEvaporationClass(true),
  maxAForFermiBreakUp(16),maxZForFermiBreakUp(8),minEForMultiFrag(1000.0*MeV) // make Multifrag. unavailable 
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

  G4FragmentVector* theResult = 0; 
  G4double exEnergy = theInitialState.GetExcitationEnergy();
  G4double A = theInitialState.GetA();
  G4int Z = theInitialState.GetZ();
  G4int Zmax = GetMaxZ();
  G4double Amax = GetMaxA();
  
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


  G4Fragment theExcitedNucleus;
  G4FragmentVector* theTempResult = 0; 

  // Check if there are excited fragments
  G4int i;
  for (i = 0; i < theResult->entries(); i++) {
    exEnergy = theResult->at(i)->GetExcitationEnergy();
    if (exEnergy > 0.0) {
      A = theResult->at(i)->GetA();
      Z = theResult->at(i)->GetZ();
      theExcitedNucleus = *(theResult->at(i));
      // try to de-excite this fragment
      if(A<GetMaxA()&&Z<GetMaxZ()&&
         exEnergy>G4NucleiPropertiesTable::GetBindingEnergy(Z,A)) {

	theTempResult = theFermiModel->BreakItUp(theExcitedNucleus);

      } else   if(exEnergy>GetMinE()*A) {


	theTempResult = theMultiFragmentation->BreakItUp(theExcitedNucleus);

      } else {

	theTempResult = theEvaporation->BreakItUp(theExcitedNucleus);

      }
      // The Nucleus has been fragmented?
      if (theTempResult->entries() > 1) {
	// If so :

	// Remove excited fragment from the result 
	delete theResult->removeAt(i--);

	// and add theTempResult elements to theResult
	while (theTempResult->entries() > 0) 
	  theResult->insert(theTempResult->removeFirst());
	  
	delete theTempResult;
      } else { // If not :
	// it doesn't matter, we Follow with the next fragment but
	// I have to make 
	theTempResult->clearAndDestroy();
	delete theTempResult;
      }
    }
  }

  //  if (theTempResult != 0 ) 
  //    {
  //      theTempResult->clearAndDestroy();
  //      delete theTempResult;
  //    }


  // Now we try to deexcite by means of PhotonEvaporation those fragments
  // which are excited.
  // In next version the Photon Evaporation has to be integrated in the main loop

  theTempResult = 0; 
  for (i = 0; i < theResult->entries(); i++) {
    if (theResult->at(i)->GetExcitationEnergy() > 0.0 && theResult->at(i)->GetA() > 1) {
      theExcitedNucleus = *(theResult->at(i));
      
      theTempResult = thePhotonEvaporation->BreakItUp(theExcitedNucleus);
      
      // Remove excited fragment from the result 
      delete theResult->removeAt(i);

      // and add theTempResult elements to theResult
      while (theTempResult->entries() > 0) 
	theResult->insert(theTempResult->removeFirst());

      theTempResult->clearAndDestroy();
      delete theTempResult;
    }
  }
  
  for (i = 0; i < theResult->entries(); i++) 
	G4LorentzVector mom(theResult->at(i)->GetMomentum());
 
    
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

  for (G4int i = 0; i < theFragmentVector->entries(); i++) {
//     theFragmentVector->at(i)->DumpInfo();
    theFragmentA = theFragmentVector->at(i)->GetA();
    theFragmentZ = theFragmentVector->at(i)->GetZ();
    theFragmentMomentum = theFragmentVector->at(i)->GetMomentum();
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
      theNew->SetFormationTime(theFragmentVector->at(i)->GetCreationTime());
      theReactionProductVector->insert(theNew);
    }
  }
  if (theFragmentVector != 0)
    { 
      theFragmentVector->clearAndDestroy();
      delete theFragmentVector;
    }
  return theReactionProductVector;
}







