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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULES:              G4NuclearDecayChannel.cc
//
// Version:             0.b.4
// Date:                14/04/00
// Author:              F Lei & P R Truscott
// Organisation:        DERA UK
// Customer:            ESA/ESTEC, NOORDWIJK
// Contract:            12115/96/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 29 February 2000, P R Truscott, DERA UK
// 0.b.3 release.
//
// 18 October 2002, F Lei
//            modified link metheds in DecayIt() to G4PhotoEvaporation() in order to
//            use the new Internal Coversion feature.      
// 13 April 2000, F Lei, DERA UK
//            Changes made are:
//            1) Use PhotonEvaporation instead of DiscreteGammaDeexcitation
//            2) verbose control
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4NuclearLevelManager.hh"
#include "G4NuclearLevelStore.hh"
#include "G4NuclearDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4IonTable.hh"

#include "G4BetaFermiFunction.hh"
#include "G4PhotonEvaporation.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4AtomicDeexcitation.hh"

const G4double G4NuclearDecayChannel:: pTolerance = 0.001;
const G4double G4NuclearDecayChannel:: levelTolerance = 2.0*keV;
//const G4bool G4NuclearDecayChannel:: FermiOn = true;

///////////////////////////////////////////////////////////////////////////////
//
//
// Constructor for one decay product (the nucleus).
//
G4NuclearDecayChannel::G4NuclearDecayChannel
                      (const G4RadioactiveDecayMode &theMode,
                       G4int Verbose,
                       const G4ParticleDefinition *theParentNucleus,
                       G4double theBR,
                       G4double theQtransition,
                       G4int A,
                       G4int Z,
                       G4double theDaughterExcitation) :
  G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<"G4NuclearDecayChannel constructor for " <<G4int(theMode) <<G4endl;}
#endif
  SetParent(theParentNucleus);
  FillParent();
  parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (1);
  FillDaughterNucleus (0, A, Z, theDaughterExcitation);
  Qtransition = theQtransition;
  halflifethreshold = 1e-6*second;
  applyICM = true;
  applyARM = true;
}
///////////////////////////////////////////////////////////////////////////////
//
//
// Constructor for a daughter nucleus and one other particle.
//
G4NuclearDecayChannel::G4NuclearDecayChannel
                      (const G4RadioactiveDecayMode &theMode,
                       G4int Verbose,
                       const G4ParticleDefinition *theParentNucleus,
                       G4double theBR,
                       G4double theQtransition,
                       G4int A,
                       G4int Z,
                       G4double theDaughterExcitation,
                       const G4String theDaughterName1) :
			G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<"G4NuclearDecayChannel constructor for " <<G4int(theMode) <<G4endl;}
#endif
  SetParent (theParentNucleus);
  FillParent();
  parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (2);
  SetDaughter(0, theDaughterName1);
  FillDaughterNucleus (1, A, Z, theDaughterExcitation);
  Qtransition = theQtransition;
  halflifethreshold = 1e-6*second;
  applyICM = true;
  applyARM = true;
}

//
// Constructor for a daughter nucleus and two other particles
//
G4NuclearDecayChannel::G4NuclearDecayChannel
                      (const G4RadioactiveDecayMode &theMode,
                       G4int Verbose,
                       const G4ParticleDefinition *theParentNucleus,
                       G4double theBR,
                       G4double /* theFFN */,
		       G4bool /* betaS */, 
		       CLHEP::RandGeneral* randBeta,
                       G4double theQtransition,
                       G4int A,
                       G4int Z,
                       G4double theDaughterExcitation,
                       const G4String theDaughterName1,
                       const G4String theDaughterName2)
 : G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<"G4NuclearDecayChannel constructor for " <<G4int(theMode) <<G4endl;}
#endif
  SetParent (theParentNucleus);
  FillParent();
  parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (3);
  SetDaughter(0, theDaughterName1);
  SetDaughter(2, theDaughterName2);
  FillDaughterNucleus(1, A, Z, theDaughterExcitation);
  RandomEnergy = randBeta;
  Qtransition = theQtransition;
  halflifethreshold = 1e-6*second;
  applyICM = true;
  applyARM = true;
}

////////////////////////////////////////////////////////////////////////////////
//
#include "G4HadTmpUtil.hh"

void G4NuclearDecayChannel::FillDaughterNucleus (G4int index, G4int A, G4int Z,
  G4double theDaughterExcitation)
{
  //
  //
  // Determine if the proposed daughter nucleus has a sensible A, Z and excitation
  // energy.
  //
  if (A<1 || Z<0 || theDaughterExcitation <0.0)
  {
    G4cerr <<"Error in G4NuclearDecayChannel::FillDaughterNucleus";
    G4cerr <<"Inappropriate values of daughter A, Z or excitation" <<G4endl;
    G4cerr <<"A = " <<A <<" and Z = " <<Z;
    G4cerr <<" Ex = " <<theDaughterExcitation*MeV  <<"MeV" <<G4endl;
    G4Exception(__FILE__, G4inttostring(__LINE__), FatalException, "G4NuclearDecayChannel::FillDaughterNucleus");
  }
  //
  //
  // Save A and Z to local variables.  Find the GROUND STATE of the daughter
  // nucleus and save this, as an ion, in the array of daughters.
  //
  daughterA = A;
  daughterZ = Z;
  if (Z == 1 && A == 1) {
    daughterNucleus = G4Proton::Definition();
  } else if (Z == 0 && A == 1) {
    daughterNucleus = G4Neutron::Definition();
  } else {
    G4IonTable *theIonTable =
      (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
    daughterNucleus = theIonTable->GetIon(daughterZ, daughterA, theDaughterExcitation*MeV);
  }
  daughterExcitation = theDaughterExcitation;
  SetDaughter(index, daughterNucleus);
}

///////////////////////////////////////////////////////////////////////////////
//
//
//
//
G4DecayProducts *G4NuclearDecayChannel::DecayIt (G4double theParentMass)
{
  //
  //
  // Load-up the details of the parent and daughter particles if they have not
  // been defined properly.
  //
  if (parent == 0) FillParent();
  if (daughters == 0) FillDaughters();
  //
  //
  // We want to ensure that the difference between the total
  // parent and daughter masses equals the energy liberated by the transition.
  //
  theParentMass = 0.0;
  for( G4int index=0; index < numberOfDaughters; index++)
    {theParentMass += daughters[index]->GetPDGMass();}
  theParentMass += Qtransition  ;
  // bug fix for beta+ decay (flei 25/09/01)
  if (decayMode == 2) theParentMass -= 2*0.511 * MeV;
  //
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1) {
    G4cout << "G4NuclearDecayChannel::DecayIt "<< G4endl;
    G4cout << "the decay mass = " << theParentMass << G4endl;
  }
#endif

  SetParentMass (theParentMass);
  
  //
  //
  // Define a product vector.
  //
  G4DecayProducts *products = 0;
  //
  //
  // Depending upon the number of daughters, select the appropriate decay
  // kinematics scheme.
  //
  switch (numberOfDaughters)
    {
    case 0:
      G4cerr << "G4NuclearDecayChannel::DecayIt ";
      G4cerr << " daughters not defined " <<G4endl;
      break;
    case 1:
      products =  OneBodyDecayIt();
      break;
    case 2:
      products =  TwoBodyDecayIt();
      break;
    case 3:
      products =  BetaDecayIt();
      break;
    default:
      G4cerr <<"Error in G4NuclearDecayChannel::DecayIt" <<G4endl;
      G4cerr <<"Number of daughters in decay = " <<numberOfDaughters <<G4endl;
      G4Exception(__FILE__, G4inttostring(__LINE__), FatalException,  "G4NuclearDecayChannel::DecayIt");
    }
  if (products == 0) {
    G4cerr << "G4NuclearDecayChannel::DecayIt ";
    G4cerr << *parent_name << " can not decay " << G4endl;
    DumpInfo();
  }
  //
  // If the decay is to an excited state of the daughter nuclide, we need
  // to apply the photo-evaporation process. This includes the IT decay mode itself.
  //
  // needed to hold the shell idex after ICM
  G4int shellIndex = -1;
  //
  if (daughterExcitation > 0.0) 
    {
      //
      // Pop the daughter nucleus off the product vector - we need to retain
      // the momentum of this particle.
      //
      dynamicDaughter = products->PopProducts();
      G4LorentzVector daughterMomentum = dynamicDaughter->Get4Momentum();
      G4ThreeVector const daughterMomentum1(static_cast<const G4LorentzVector> (daughterMomentum));
      //
      //
      // Now define a G4Fragment with the correct A, Z and excitation, and declare and
      // initialise a G4PhotonEvaporation object.
      //    
      G4Fragment nucleus(daughterA, daughterZ, daughterMomentum);
      G4PhotonEvaporation* deexcitation = new G4PhotonEvaporation;
      deexcitation->SetVerboseLevel(GetVerboseLevel());
      // switch on/off internal electron conversion
      deexcitation->SetICM(applyICM);
      // set the maximum life-time for a level that will be treated. Level with life-time longer than this
      // will be outputed as meta-stable isotope
      //
      deexcitation->SetMaxHalfLife(halflifethreshold);
      // but in IT mode, we need to force the transition 
      if (decayMode == 0) {
	deexcitation->RDMForced(true);
      } else {
	deexcitation->RDMForced(false);
      }
      //
      // Now apply the photo-evaporation
      // Use BreakUp() so limit to one transition at a time, if ICM is requested
      // this change is realted to bug#1001  (F.Lei 07/05/2010)
      G4FragmentVector* gammas = 0;	
      if (applyICM) {
	gammas = deexcitation->BreakUp(nucleus);	
      } else {
	gammas = deexcitation->BreakItUp(nucleus);
      }
      // the returned G4FragmentVector contains the residual nuclide
      // as its last entry.
      G4int nGammas=gammas->size()-1;
      //
      // Go through each gamma/e- and add it to the decay product.  The angular distribution
      // of the gammas is isotropic, and the residual nucleus is assumed not to have suffered
      // any recoil as a result of this de-excitation.
      //
      for (G4int ig=0; ig<nGammas; ig++)
	{
	  G4DynamicParticle *theGammaRay = new
	    G4DynamicParticle (gammas->operator[](ig)->GetParticleDefinition(),
			       gammas->operator[](ig)->GetMomentum());
	  theGammaRay -> SetProperTime(gammas->operator[](ig)->GetCreationTime());
	  products->PushProducts (theGammaRay);
	}
      //
      // now the nucleus
      G4double finalDaughterExcitation = gammas->operator[](nGammas)->GetExcitationEnergy();
      // f.lei (03/01/03) this is needed to fix the crach in test18 
      if (finalDaughterExcitation <= 1.0*keV) finalDaughterExcitation = 0 ;
      
      // f.lei (07/03/05) added the delete to fix bug#711
      if (dynamicDaughter) delete dynamicDaughter;
      
      G4IonTable *theIonTable =  (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());      
      dynamicDaughter = new G4DynamicParticle
	(theIonTable->GetIon(daughterZ,daughterA,finalDaughterExcitation),
	 daughterMomentum1);
      products->PushProducts (dynamicDaughter); 
      
      // retrive the ICM shell index
      shellIndex = deexcitation->GetVacantShellNumber();
      
      //
      // Delete/reset variables associated with the gammas.
      //
      while (!gammas->empty()) {
	delete *(gammas->end()-1);
	gammas->pop_back();
      }
      //    gammas->clearAndDestroy();
      delete gammas;
      delete deexcitation;
    }
  //
  // now we have to take care of the EC product which have to go through the ARM
  // 
  G4int eShell = -1;
  if (decayMode == 3 || decayMode == 4 || decayMode == 5) {
    switch (decayMode)
      {
      case KshellEC:
	//
	{
	  eShell = 0; // --> 0 from 1 (f.lei 30/4/2008)
	}
	break;
      case LshellEC:
	//
	{
	  eShell = G4int(G4UniformRand()*3)+1;
	}
	break;
      case MshellEC:
	//
	{
	// limit the shell index to 6 as specified by the ARM (F.Lei 06/05/2010)
	// eShell = G4int(G4UniformRand()*5)+4;
	  eShell = G4int(G4UniformRand()*3)+4;
	}
	break;
      case ERROR:
      default:
        G4Exception("G4NuclearDecayChannel::DecayIt()", "601",
                    FatalException, "Error in decay mode selection");
      }
  }
  // Also have to deal with the IT case where ICM may have been applied
  //
  if (decayMode == 0) {
    eShell = shellIndex;
  }
  //
  // now apply ARM if it is requested and there is a vaccancy
  //
  if (applyARM && eShell != -1) {
    G4int aZ = daughterZ;
    if (aZ > 5 && aZ < 100) {  // only applies to 5< Z <100 
      // Retrieve the corresponding identifier and binding energy of the selected shell
      const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();

      //check that eShell is smaller than the max number of shells
      //this to avoid a warning from transitionManager(otherwise the output is the same)
      //Correction to Bug 1662. L Desorgher (04/02/2011)
      if (eShell >= transitionManager->NumberOfShells(aZ)){
    	  eShell=transitionManager->NumberOfShells(aZ)-1;
      }

      const G4AtomicShell* shell = transitionManager->Shell(aZ, eShell);
      G4double bindingEnergy = shell->BindingEnergy();
      G4int shellId = shell->ShellId();

      G4AtomicDeexcitation* atomDeex = new G4AtomicDeexcitation();  
      //the default is no Auger electron generation. 
      // Switch it on/off here! 
      atomDeex->ActivateAugerElectronProduction(true);
      std::vector<G4DynamicParticle*>* armProducts = atomDeex->GenerateParticles(aZ,shellId);

      // pop up the daughter before insertion; 
      // f.lei (30/04/2008) check if the total kinetic energy is less than 
      // the shell binding energy; if true add the difference to the daughter to conserve the energy 
      dynamicDaughter = products->PopProducts();
      G4double tARMEnergy = 0.0; 
      for (size_t i = 0;  i < armProducts->size(); i++) {
	products->PushProducts ((*armProducts)[i]);
	tARMEnergy += (*armProducts)[i]->GetKineticEnergy();
      }
      if ((bindingEnergy - tARMEnergy) > 0.1*keV){
	G4double dEnergy = dynamicDaughter->GetKineticEnergy() + (bindingEnergy - tARMEnergy);
	dynamicDaughter->SetKineticEnergy(dEnergy);
      }
      products->PushProducts(dynamicDaughter); 

#ifdef G4VERBOSE
      if (GetVerboseLevel()>0)
	{
      	  G4cout <<"G4NuclearDecayChannel::Selected shell number for ARM =  " <<shellId <<G4endl;
	  G4cout <<"G4NuclearDecayChannel::ARM products =  " <<armProducts->size()<<G4endl;
	  G4cout <<"                 The binding energy =  " << bindingEnergy << G4endl;              
	  G4cout <<"  Total ARM particle kinetic energy =  " << tARMEnergy << G4endl;              
	}
#endif   

      delete armProducts;
      delete atomDeex;
    }
  }

  return products;
}


G4DecayProducts* G4NuclearDecayChannel::BetaDecayIt()
{
  if (GetVerboseLevel()>1) G4cout << "G4Decay::BetaDecayIt()"<<G4endl;

  //daughters'mass
  G4double daughtermass[3];
  G4double sumofdaughtermass = 0.0;
  G4double pmass = GetParentMass();
  for (G4int index=0; index<3; index++)
    {
     daughtermass[index] = daughters[index]->GetPDGMass();
     sumofdaughtermass += daughtermass[index];
    }

  //create parent G4DynamicParticle at rest
  G4ParticleMomentum dummy;
  G4DynamicParticle * parentparticle = new G4DynamicParticle( parent, dummy, 0.0);

  //create G4Decayproducts
  G4DecayProducts *products = new G4DecayProducts(*parentparticle);
  delete parentparticle;

  G4double Q = pmass - sumofdaughtermass;  

  // faster method as suggested by Dirk Kruecker of FZ-Julich
  G4double daughtermomentum[3];
  G4double daughterenergy[3];
  // Use the histogram distribution to generate the beta energy
  daughterenergy[0] = RandomEnergy->shoot() * Q;
  daughtermomentum[0] = std::sqrt(daughterenergy[0]*daughterenergy[0] +
                        2.0*daughterenergy[0] * daughtermass[0]);

  // neutrino energy distribution is flat within the kinematical limits
  G4double rd = 2*G4UniformRand()-1;
  // limits
  G4double Mme=pmass-daughtermass[0];
  G4double K=0.5-daughtermass[1]*daughtermass[1]/(2*Mme*Mme-4*pmass*daughterenergy[0]);
	  
  daughterenergy[2]=K*(Mme-daughterenergy[0]+rd*daughtermomentum[0]);
  daughtermomentum[2] = daughterenergy[2] ; 
	  
  // the recoil nucleus
  daughterenergy[1] = Q-daughterenergy[0]-daughterenergy[2];
  G4double recoilmomentumsquared = daughterenergy[1]*daughterenergy[1] +
                             2.0*daughterenergy[1] * daughtermass[1];
  if (recoilmomentumsquared < 0.0) recoilmomentumsquared = 0.0;
  daughtermomentum[1] = std::sqrt(recoilmomentumsquared);
  
  // output message
  if (GetVerboseLevel()>1) {
    G4cout <<"     daughter 0:" <<daughtermomentum[0]/GeV <<"[GeV/c]" <<G4endl;
    G4cout <<"     daughter 1:" <<daughtermomentum[1]/GeV <<"[GeV/c]" <<G4endl;
    G4cout <<"     daughter 2:" <<daughtermomentum[2]/GeV <<"[GeV/c]" <<G4endl;
  }
  //create daughter G4DynamicParticle
  G4double costheta, sintheta, phi, sinphi, cosphi;
  G4double costhetan, sinthetan, phin, sinphin, cosphin;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
  G4ParticleMomentum direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle
      = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);
    
  costhetan = (daughtermomentum[1]*daughtermomentum[1]-
               daughtermomentum[2]*daughtermomentum[2]-
               daughtermomentum[0]*daughtermomentum[0])/
        (2.0*daughtermomentum[2]*daughtermomentum[0]);
  // added the following test to avoid rounding erros. A problem
  // reported bye Ben Morgan of Uni.Warwick
  if (costhetan > 1.) costhetan = 1.;
  if (costhetan < -1.) costhetan = -1.;
  sinthetan = std::sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = twopi*G4UniformRand()*rad;
  sinphin = std::sin(phin);
  cosphin = std::cos(phin);
  G4ParticleMomentum direction2;
  direction2.setX(sinthetan*cosphin*costheta*cosphi - 
                  sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi);
  direction2.setY(sinthetan*cosphin*costheta*sinphi +
                  sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi);
  direction2.setZ(-sinthetan*cosphin*sintheta + costhetan*costheta);
  daughterparticle = new G4DynamicParticle(daughters[2],
                          direction2*(daughtermomentum[2]/direction2.mag()));
  products->PushProducts(daughterparticle);
    
  daughterparticle =
    new G4DynamicParticle(daughters[1],
                         (direction0*daughtermomentum[0] +
			  direction2*(daughtermomentum[2]/direction2.mag()))*(-1.0));
  products->PushProducts(daughterparticle);
  
  if (GetVerboseLevel()>1) {
    G4cout << "G4NuclearDecayChannel::BetaDecayIt ";
    G4cout << "  create decay products in rest frame " <<G4endl;
    products->DumpInfo();
  }
  return products;
}
