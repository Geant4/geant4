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
// 17 October 2011, L. Desorgher
//			  -Allow the atomic relaxation after de-excitation of exited 
//                          nuclei even for beta and alpha
//			    decay. Bug found and solution proposed by Ko Abe.
//			  -Set halflifethreshold by default to a negative value
//
// 20 November 2011, V.Ivanchenko
//                        - Migration to new design of atomic deexcitation
//
///////////////////////////////////////////////////////////////////////////////

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
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

#include "G4VAtomDeexcitation.hh"
#include "G4AtomicShells.hh"
#include "G4LossTableManager.hh"

//Model const parameters
const G4double G4NuclearDecayChannel:: pTolerance = 0.001;
const G4double G4NuclearDecayChannel:: levelTolerance = 2.0*keV;
//const G4bool G4NuclearDecayChannel:: FermiOn = true;

//This is a kind of "cache"
G4ThreadLocal G4DynamicParticle* G4NuclearDecayChannel::dynamicDaughter = 0;

// Constructor for one decay product (the nucleus)
G4NuclearDecayChannel::
G4NuclearDecayChannel(const G4RadioactiveDecayMode& theMode,
                      G4int Verbose,
                      const G4ParticleDefinition* theParentNucleus,
                      const G4double theBR,
                      const G4double theQtransition,
                      const G4int A, const G4int Z,
                      const G4double theDaughterExcitation)
 :G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode),
  Qtransition(theQtransition), RandomEnergy(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "G4NuclearDecayChannel constructor for " << G4int(theMode)
           << G4endl;
  }
#endif
  SetParent(theParentNucleus);
  FillParent();
  G4MT_parent_mass = theParentNucleus->GetPDGMass();
  SetBR(theBR);
  SetNumberOfDaughters (1);
  FillDaughterNucleus(0, A, Z, theDaughterExcitation);
  halflifethreshold = nanosecond;
  applyICM = true;
  applyARM = true;
  FillDaughters();
}

// Constructor for a daughter nucleus and one other particle.
//
G4NuclearDecayChannel::
G4NuclearDecayChannel(const G4RadioactiveDecayMode& theMode,
                      G4int Verbose,
                      const G4ParticleDefinition* theParentNucleus,
                      const G4double theBR,
                      const G4double theQtransition,
                      const G4int A, const G4int Z,
                      const G4double theDaughterExcitation,
                      const G4String theDaughterName1)
 :G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode),
  Qtransition(theQtransition), RandomEnergy(0)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "G4NuclearDecayChannel constructor for " << G4int(theMode)
           << G4endl;
  }
#endif
  SetParent(theParentNucleus);
  FillParent();
  G4MT_parent_mass = theParentNucleus->GetPDGMass();
  SetBR(theBR);
  SetNumberOfDaughters (2);
  SetDaughter(0, theDaughterName1);
  FillDaughterNucleus(1, A, Z, theDaughterExcitation);
  halflifethreshold = nanosecond;
  applyICM = true;
  applyARM = true;
  FillDaughters();
}

// Constructor for a daughter nucleus and two other particles
//
G4NuclearDecayChannel::
G4NuclearDecayChannel(const G4RadioactiveDecayMode &theMode,
                      G4int Verbose,
                      const G4ParticleDefinition *theParentNucleus,
                      const G4double theBR,
                      G4double /* theFFN */,
		      G4bool /* betaS */, 
		      G4RandGeneral* randBeta,
                      const G4double theQtransition,
                      const G4int A, const G4int Z,
                      const G4double theDaughterExcitation,
                      const G4String theDaughterName1,
                      const G4String theDaughterName2)
 :G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode),
  Qtransition(theQtransition)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel() > 1) {
    G4cout << "G4NuclearDecayChannel constructor for " << G4int(theMode)
           << G4endl;
  }
#endif
  SetParent(theParentNucleus);
  FillParent();
  G4MT_parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (3);
  SetDaughter(0, theDaughterName1);
  SetDaughter(2, theDaughterName2);
  FillDaughterNucleus(1, A, Z, theDaughterExcitation);
  RandomEnergy = randBeta;

  halflifethreshold = nanosecond;
  applyICM = true;
  applyARM = true;
  FillDaughters();
}

G4NuclearDecayChannel::~G4NuclearDecayChannel()
{} 

void G4NuclearDecayChannel::FillDaughterNucleus(G4int index, G4int A, G4int Z,
                                                const G4double theDaughterExcitation)
{
  // Determine if the proposed daughter nucleus has a sensible A, Z and
  // excitation energy.
  if (A < 1 || Z < 0 || theDaughterExcitation < 0.0) {
    G4ExceptionDescription ed;
    ed << "Inappropriate values of daughter A, Z or excitation: "
       << A << " , " << Z << " , " << theDaughterExcitation*MeV << " MeV "
       << G4endl;
    G4Exception("G4NuclearDecayChannel::FillDaughterNucleus()", "HAD_RDM_006",
                FatalException, ed);
  }

  // Save A and Z to local variables.  Find the GROUND STATE of the daughter
  // nucleus and save this, as an ion, in the array of daughters.
  daughterA = A;
  daughterZ = Z;
  if (Z == 1 && A == 1) {
    daughterNucleus = G4Proton::Definition();
  } else if (Z == 0 && A == 1) {
    daughterNucleus = G4Neutron::Definition();
  } else {
    G4IonTable *theIonTable =
      (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
    // GetIon with only Z and A arguments returns ground state         
    daughterNucleus = theIonTable->GetIon(daughterZ, daughterA);
  }

  daughterExcitation = theDaughterExcitation;
  SetDaughter(index, daughterNucleus);
}

G4DecayProducts* G4NuclearDecayChannel::DecayIt(G4double)
{
  G4double deltaM;
  if (decayMode == 1) {                         // beta- decay
    deltaM = CLHEP::electron_mass_c2;
  } else if (decayMode == 2) {                  // beta+ decay
    deltaM = 2.*CLHEP::electron_mass_c2;
  } else if (decayMode < 6 && decayMode > 2) {  // EC
    deltaM = -CLHEP::electron_mass_c2;
  } else {                                      // all others
    deltaM = 0.0;
  }

  // Mass available for decay in rest frame of parent after correcting for
  // the appropriate number of electron masses and reserving the daughter
  // excitation energy to be applied later 
  G4double massOfParent = G4MT_parent->GetPDGMass();  // PDG mass includes excitation energy
  SetParentMass(massOfParent - deltaM - daughterExcitation);

  // Define a product vector.
  G4DecayProducts* products = 0;

  // Depending upon the number of daughters, select the appropriate decay
  // kinematics scheme.
  switch (numberOfDaughters) {
    case 0:
      {
      G4ExceptionDescription ed;
      ed << " No daughters defined " << G4endl;
      G4Exception("G4NuclearDecayChannel::DecayIt()", "HAD_RDM_005",
                  JustWarning, ed);
      }
      break;
    case 1:
      products = OneBodyDecayIt();
      break;
    case 2:
      products = TwoBodyDecayIt();
      break;
    case 3:
      products = BetaDecayIt();
      break;
    default:
    {
      G4ExceptionDescription ed;
      ed << " More than 3 daughters in decay: N = " << numberOfDaughters
         << G4endl;
      G4Exception("G4NuclearDecayChannel::DecayIt()", "HAD_RDM_007",
                  FatalException, ed);
    }
  }

  if (products == 0) {
    G4ExceptionDescription ed;
    ed << " Parent nucleus " << *parent_name << " was not decayed " << G4endl;
    G4Exception("G4NuclearDecayChannel::DecayIt()", "HAD_RDM_008",
                JustWarning, ed);
    DumpInfo();
  } else {
    // If the decay is to an excited state of the daughter nuclide, the photon
    // evaporation process must be applied.

    // Need to hold the shell idex after ICM
    G4int shellIndex = -1;
    if (daughterExcitation > 0.0) {
      // Pop the daughter nucleus off the product vector to get its 4-momentum
      dynamicDaughter = products->PopProducts();
      G4LorentzVector daughterMomentum = dynamicDaughter->Get4Momentum();
      if (dynamicDaughter) delete dynamicDaughter;

      // Using daughter nucleus, set up a G4Fragment for photon evaporatation
      if (decayMode == 0) {
        G4double exe = ((const G4Ions*)(G4MT_parent))->GetExcitationEnergy();
        daughterMomentum.setE(daughterMomentum.e() + exe);
      }
      G4Fragment nucleus(daughterA, daughterZ, daughterMomentum);
      G4PhotonEvaporation* deexcitation = new G4PhotonEvaporation;
      deexcitation->SetVerboseLevel(GetVerboseLevel());

      // switch on/off internal electron conversion
      deexcitation->SetICM(applyICM);

      // In IT mode, we need to force the transition 
      if (decayMode == 0) {
        deexcitation->RDMForced(true);
        // at this point, de-excitation will occur even if IT state is long-lived (>1ns)
        // Why does it need to be forced? 
      } else {
        // Not forced, but decay will still happen if lifetime < 1 ns,
        // otherwise no gamma decay is performed
        deexcitation->RDMForced(false);
      }

      //////////////////////////////////////////////////////////////////////////
      //                                                                      //
      //  Apply photon evaporation if Isomeric Transition is indicated.       //
      //  Use G4PhotonEvaporation::BreakUp() which does only one transition.  //
      //  This allows IC to be done.                                          //
      //                                                                      //
      //////////////////////////////////////////////////////////////////////////

      G4IonTable* theIonTable =
              (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
      G4ParticleDefinition* daughterIon = 0;

      if (decayMode != 0) {
        daughterIon = theIonTable->GetIon(daughterZ, daughterA, daughterExcitation);
      } else {
        // The fragment vector from photon evaporation contains the list of
        // evaporated gammas, some of which may have been replaced by conversion
        // electrons.  The last element is the residual nucleus.
        // Note: try photoEvapProducts as a name instead of gammas
        G4FragmentVector* gammas = deexcitation->BreakUp(nucleus);
        G4int nFrags = G4int(gammas->size());
        G4double eOrGammaEnergy = 0.0;

        if (nFrags < 1) {
          G4ExceptionDescription ed;
          ed << nFrags << " No fragments produced by photon evaporation. " << G4endl;
          G4Exception("G4NuclearDecayChannel::DecayIt()","HAD_RDM_012",
                      FatalException, ed);
        } else if (nFrags > 1) {
          // Add gamma/e- to the decay product. The angular distribution of this 
          // particle is assumed to be isotropic
          G4Fragment* eOrGamma;
          G4DynamicParticle* eOrGammaDyn;
          for (G4int i = 0; i < nFrags - 1; i++) {
            eOrGamma = gammas->operator[](i);
            eOrGammaDyn = new G4DynamicParticle(eOrGamma->GetParticleDefinition(),
                                                eOrGamma->GetMomentum() );
            eOrGammaDyn->SetProperTime(eOrGamma->GetCreationTime() );
            products->PushProducts(eOrGammaDyn);
            eOrGammaEnergy += eOrGamma->GetMomentum().e();
          }
        } 

        G4double finalDaughterExcitation =
                            gammas->operator[](nFrags-1)->GetExcitationEnergy();
        if (finalDaughterExcitation <= 1.0*keV) finalDaughterExcitation = 0;

        // Get new ion with excitation energy reduced by emitted gamma energy
        daughterIon =
             theIonTable->GetIon(daughterZ, daughterA, finalDaughterExcitation);
        daughterMomentum.setE(daughterMomentum.e() - eOrGammaEnergy);

        // Delete/reset variables associated with the gammas.
        while (!gammas->empty() ) {
          delete *(gammas->end()-1);
          gammas->pop_back();
        }
        delete gammas;
      } // end if decayMode == 0
 
      G4ThreeVector const daughterMomentum1(static_cast<const G4LorentzVector> (daughterMomentum));
      dynamicDaughter = new G4DynamicParticle(daughterIon, daughterMomentum1);  
      products->PushProducts(dynamicDaughter);

      // retrieve the ICM shell index
      shellIndex = deexcitation->GetVacantShellNumber();
      
      delete deexcitation;
    } // if daughter excitation > 0

    // Now take care of the EC products which have to go through the ARM
    G4int eShell = -1;
    if (decayMode == 3 || decayMode == 4 || decayMode == 5) {
      switch (decayMode)
        {
        case KshellEC:
	  {
            eShell = 0; // --> 0 from 1 (f.lei 30/4/2008)
          }
          break;
        case LshellEC:
          {
            eShell = G4int(G4UniformRand()*3)+1;
          }
          break;
        case MshellEC:
          {
            // limit the shell index to 6 as specified by the ARM (F.Lei 06/05/2010)
            // eShell = G4int(G4UniformRand()*5)+4;
            eShell = G4int(G4UniformRand()*3)+4;
          }
          break;
        case RDM_ERROR:
        default:
        G4Exception("G4NuclearDecayChannel::DecayIt()", "HAD_RDM_009",
                    FatalException, "Incorrect decay mode selection");
        }
    } else {
      // For other cases eShell comes from shellIndex resulting from  the photo decay
      // modeled by G4PhotonEvaporation* de-excitation (see above)
      eShell = shellIndex;
    }

    // now apply ARM if it is requested and there is a vaccancy
    if (applyARM && eShell != -1) {
      G4int aZ = daughterZ;

      // V.Ivanchenko migration to new interface to atomic deexcitation
      // no check on index of G4MaterialCutsCouple, simplified 
      // check on secondary energy Esec < 0.1 keV
      G4VAtomDeexcitation* atomDeex =
               G4LossTableManager::Instance()->AtomDeexcitation();
      if (atomDeex) {
        if(atomDeex->IsFluoActive() && aZ > 5 && aZ < 100) {  // only applies to 5< Z <100 
          if (eShell >= G4AtomicShells::GetNumberOfShells(aZ)){
            eShell = G4AtomicShells::GetNumberOfShells(aZ)-1;
          }
          G4AtomicShellEnumerator as = G4AtomicShellEnumerator(eShell);
          const G4AtomicShell* shell = atomDeex->GetAtomicShell(aZ, as);    
          std::vector<G4DynamicParticle*> armProducts;
          const G4double deexLimit = 0.1*keV;
          atomDeex->GenerateParticles(&armProducts, shell, aZ, deexLimit, deexLimit);
          size_t narm = armProducts.size();
          if (narm > 0) {
	    // L.Desorgher: need to initialize dynamicDaughter in some decay
            // cases (for example Hg194)
	    dynamicDaughter = products->PopProducts();
	    G4ThreeVector bst = dynamicDaughter->Get4Momentum().boostVector();
            for (size_t i = 0; i<narm; ++i) {
              G4DynamicParticle* dp = armProducts[i];
              G4LorentzVector lv = dp->Get4Momentum().boost(bst);
              dp->Set4Momentum(lv);
              products->PushProducts(dp);
            }
            products->PushProducts(dynamicDaughter);
          }
        }
      }
    }
  } // Parent nucleus decayed
  /*
    if (atomDeex && aZ > 5 && aZ < 100) {  // only applies to 5< Z <100 
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
  */

  return products;
}

G4DecayProducts* G4NuclearDecayChannel::BetaDecayIt()
{
  G4double pmass = GetParentMass();

  G4double daughtermass[3];
  for (G4int index = 0; index < 3; index++) {
    daughtermass[index] = G4MT_daughters[index]->GetPDGMass();
  }

  // Add excitation energy so daughter can be decayed later
  daughtermass[1] += daughterExcitation;

  // Create parent G4DynamicParticle at rest and create products
  G4ParticleMomentum dummy;
  G4DynamicParticle parentParticle(G4MT_parent, dummy, 0.0);
  G4DecayProducts* products = new G4DecayProducts(parentParticle);

  // faster method as suggested by Dirk Kruecker of FZ-Julich
  G4double daughtermomentum[3];
  G4double daughterenergy[3];
  // Use the histogram distribution to generate the beta energy
  // 0 = electron, 1 = daughter, 2 = neutrino
  daughterenergy[0] = Qtransition*RandomEnergy->shoot(G4Random::getTheEngine());
  daughtermomentum[0] = std::sqrt(daughterenergy[0]*(daughterenergy[0] + 2.*daughtermass[0]) );

  // neutrino energy distribution is flat within the kinematical limits
  G4double rd = 2.*G4UniformRand() - 1.;
  // limits
  G4double Mme = daughtermass[1] + Qtransition;
  G4double K = 0.5 - daughtermass[1]*daughtermass[1]/(2*Mme*Mme-4*pmass*daughterenergy[0]);
	  
  daughterenergy[2] = K * (Mme - daughterenergy[0] + rd*daughtermomentum[0]);
  daughtermomentum[2] = daughterenergy[2];

  // the recoil nucleus
  daughterenergy[1] = Qtransition - daughterenergy[0] - daughterenergy[2];
  G4double recoilmomentumsquared =
               daughterenergy[1]*(daughterenergy[1] + 2.0*daughtermass[1]);
  if (recoilmomentumsquared < 0.0) recoilmomentumsquared = 0.0;
  daughtermomentum[1] = std::sqrt(recoilmomentumsquared);
  
  // output message
  if (GetVerboseLevel()>1) {
    G4cout << " G4NuclearDecayChannel::BetaDecayIt() " << G4endl;
    G4cout <<"     e- momentum: " <<daughtermomentum[0]/GeV <<" [GeV/c]" <<G4endl;
    G4cout <<"     daughter momentum: " <<daughtermomentum[1]/GeV <<" [GeV/c]" <<G4endl;
    G4cout <<"     nu momentum: " <<daughtermomentum[2]/GeV <<" [GeV/c]" <<G4endl;
    G4cout <<"     e- energy: " << daughtermass[0] + daughterenergy[0] << G4endl;
    G4cout <<"     daughter energy: " << daughtermass[1] + daughterenergy[1] << G4endl;
    G4cout <<"     nu energy: " << daughtermass[2] + daughterenergy[2] << G4endl;
    G4cout <<"     total of daughter energies: " << daughtermass[0] + daughtermass[1] +
               daughtermass[2] + daughterenergy[0] + daughterenergy[1] + daughterenergy[2] 
            << G4endl; 
  }
  //create daughter G4DynamicParticle
  G4double costheta, sintheta, phi, sinphi, cosphi;
  G4double costhetan, sinthetan, phin, sinphin, cosphin;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = std::sqrt((1.0-costheta)*(1.0+costheta));
  phi  = twopi*G4UniformRand()*rad;
  sinphi = std::sin(phi);
  cosphi = std::cos(phi);
// electron chosen isotropically
  G4ParticleMomentum direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle
      = new G4DynamicParticle( G4MT_daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);
  // cos of angle between electron and neutrino
  costhetan = (daughtermomentum[1]*daughtermomentum[1]-
               daughtermomentum[2]*daughtermomentum[2]-
               daughtermomentum[0]*daughtermomentum[0])/
        (2.0*daughtermomentum[2]*daughtermomentum[0]);

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
  daughterparticle = new G4DynamicParticle(G4MT_daughters[2],
                          direction2*(daughtermomentum[2]/direction2.mag()));
  products->PushProducts(daughterparticle);
  // daughter nucleus p = - (p_e + p_nu )
  daughterparticle =
    new G4DynamicParticle(G4MT_daughters[1],
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
