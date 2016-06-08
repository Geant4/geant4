// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              G4NuclearDecayChannel.cc
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
// 13 April 2000, F Lei, DERA UK
//            Changes made are:
//            1) Use PhotonEvaporation instead of DiscreteGammaDeexcitation
//            2) verbose control
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
///////////////////////////////////////////////////////////////////////////////
//
#include "G4NuclearDecayChannel.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4ParticleChangeForRadDecay.hh"
#include "G4IonTable.hh"

#include "G4BetaFermiFunction.hh"
#include "G4DiscreteGammaDeexcitation.hh"
#include "G4PhotonEvaporation.hh"
#include "G4VGammaDeexcitation.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"

const G4double G4NuclearDecayChannel:: pTolerance = 0.001;
const G4double G4NuclearDecayChannel:: levelTolerance = 2.0*keV;
const G4bool G4NuclearDecayChannel:: FermiOn = true;

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
    {G4cout <<"G4NuclearDecayChannel constructor for " <<theMode <<G4endl;}
#endif
  SetParent(theParentNucleus);
  FillParent();
  parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (1);
  FillDaughterNucleus (0, A, Z, theDaughterExcitation);
  Qtransition = theQtransition;
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
    {G4cout <<"G4NuclearDecayChannel constructor for " <<theMode <<G4endl;}
#endif
  SetParent (theParentNucleus);
  FillParent();
  parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (2);
  SetDaughter(0, theDaughterName1);
  FillDaughterNucleus (1, A, Z, theDaughterExcitation);
  Qtransition = theQtransition;
}
///////////////////////////////////////////////////////////////////////////////
//
//
// Constructor for a daughter nucleus and two other particles.
//
G4NuclearDecayChannel::G4NuclearDecayChannel
                      (const G4RadioactiveDecayMode &theMode,
                       G4int Verbose,
                       const G4ParticleDefinition *theParentNucleus,
                       G4double theBR,
                       G4double theFFN,
                       G4double theQtransition,
                       G4int A,
                       G4int Z,
                       G4double theDaughterExcitation,
                       const G4String theDaughterName1,
                       const G4String theDaughterName2) :
  G4GeneralPhaseSpaceDecay(Verbose), decayMode(theMode)
{
#ifdef G4VERBOSE
  if (GetVerboseLevel()>1)
    {G4cout <<"G4NuclearDecayChannel constructor for " <<theMode <<G4endl;}
#endif
  SetParent (theParentNucleus);
  FillParent();
  parent_mass = theParentNucleus->GetPDGMass();
  SetBR (theBR);
  SetNumberOfDaughters (3);
  SetDaughter(0, theDaughterName1);
  SetDaughter(2, theDaughterName2);
  FillDaughterNucleus(1, A, Z, theDaughterExcitation);
  Qtransition = theQtransition;
  FermiFN = theFFN;
}

////////////////////////////////////////////////////////////////////////////////
//
//
//
//
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
    G4Exception("G4NuclearDecayChannel::FillDaughterNucleus");
  }
  //
  //
  // Save A and Z to local variables.  Find the GROUND STATE of the daughter
  // nucleus and save this, as an ion, in the array of daughters.
  //
  daughterA = A;
  daughterZ = Z;
  G4IonTable *theIonTable = (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
  //  daughterNucleus = theIonTable->GetIon(daughterZ, daughterA, 0.0*keV);
  //
  //
  // Determine the excitation state corresponds to an actual level in the
  // photo-evaporation data.  Flag an error if the difference is too large.
  //
  if (theDaughterExcitation > 0.0) {
    G4NuclearLevelManager levelManager = G4NuclearLevelManager
      (daughterZ, daughterA);
    if ( levelManager.NumberOfLevels() ) {
      const G4NuclearLevel* level = levelManager.NearestLevel (theDaughterExcitation);

      daughterExcitation = level->Energy();

      if (abs(daughterExcitation-theDaughterExcitation)>levelTolerance){
#ifdef G4VERBOSE
	if (GetVerboseLevel()>1){
	  G4cout <<"In G4NuclearDecayChannel::FillDaughterNucleus" <<G4endl;
	  G4cout <<"Difference in daughter excitation and G4NuclearLevelManager data ";
	  G4cout <<"exceeds tolerance" <<G4endl;
	  G4cout <<"Level requested = " <<theDaughterExcitation*MeV <<" MeV" <<G4endl;
	  G4cout <<"Level found     = " <<daughterExcitation*MeV <<" MeV" <<G4endl;
	  G4cout << " -- The requested energy level will be used!-- "<< G4endl;
	}
#endif
	daughterExcitation = theDaughterExcitation;
      }
      // Level hafe life is in ns and I want to set the gate as 1 micros
      // also we have to force the IT case in all conditions     	
      if (level->HalfLife() <= 1000. || index ==  0) {
	daughterNucleus = theIonTable->GetIon(daughterZ, daughterA, 0.0*keV);
      }
      else{
	daughterNucleus = theIonTable->GetIon(daughterZ, daughterA, daughterExcitation);
	daughterExcitation = 0.0;
      }     
    }
    else{
#ifdef G4VERBOSE
      if (GetVerboseLevel()>0){
	G4cout << "Error in G4NuclearDecayChannel::FillDaughterNucleus" <<G4endl;
	G4cout << "PhotonEvaporation data is not available " <<G4endl;
	G4cout << "RDM could crash during Photo De-excitaion "<< G4endl;
      } 
#endif
      daughterNucleus = theIonTable->GetIon(daughterZ, daughterA, 0.0*keV);
      daughterExcitation = theDaughterExcitation;
    }
  }
  else {
    daughterExcitation = 0.0;
    daughterNucleus = theIonTable->GetIon(daughterZ, daughterA, 0.0*keV);
  }
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
  if (parent == NULL) FillParent();
  if (daughters == NULL) FillDaughters();
  //
  //
  // THIS IS A CHEAT!  We want to ensure that the difference between the total
  // parent and daughter masses equals the energy liberated by the transition.
  //
  theParentMass = 0.0;
  for( G4int index=0; index < numberOfDaughters; index++)
    {theParentMass += daughters[index]->GetPDGMass();}
  theParentMass += Qtransition  ;

  if (GetVerboseLevel()>1) {
    G4cout << "G4NuclearDecayChannel::DecayIt ";
    G4cout << "the decay mass = " << theParentMass << G4endl;
  }
  
  SetParentMass (theParentMass);
  
  //
  //
  // Define a product vector.
  //
  G4DecayProducts *products = NULL;
  //
  //
  // Depending upon the number of daughters, select the appropriate decay
  // kinematics scheme.
  //
  switch (numberOfDaughters)
    {
    case 0:
      if (GetVerboseLevel()>0)
      {
        G4cout << "G4NuclearDecayChannel::DecayIt ";
        G4cout << " daughters not defined " <<G4endl;
      }
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
      G4Exception ("G4NuclearDecayChannel::DecayIt");
  }
  if ((products == NULL) && (GetVerboseLevel()>0)) {
    G4cerr << "G4NuclearDecayChannel::DecayIt ";
    G4cerr << *parent_name << " can not decay " << G4endl;
    DumpInfo();
  }
  //
  //
  // If the decay is to an excited state of the daughter nuclide, we need
  // to apply the photo-evaporation process.
  //
  if (daughterExcitation > 0.0)
  {
  //
  //
  // Pop the daughter nucleus off the product vector - we need to retain
  // the momentum of this particle.
  //
    dynamicDaughter = products->PopProducts();
    G4LorentzVector daughterMomentum = dynamicDaughter->Get4Momentum();
      G4ThreeVector daughterMomentum1 = dynamicDaughter->Get4Momentum();
  //
  //
  // Now define a G4Fragment with the correct A, Z and excitation, and declare and
  // initialise a G4DiscreteGammaDeexcitation object.
  //

    //  daughterMomentum.setT(daughterMomentum.t()+G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass( daughterZ, daughterA )+daughterExcitation);

    daughterMomentum.setT(daughterMomentum.t()+daughterExcitation);
    G4Fragment nucleus(daughterA, daughterZ, daughterMomentum);
    //G4LorentzVector p4(0.,0.,0.,G4NucleiProperties::GetNuclearMass(daughterA,daughterZ)
    //	       +daughterExcitation);
    //G4Fragment nucleus(daughterA, daughterZ, p4);
    //    nucleus.SetExcitationEnergy(daughterExcitation);
   

    // G4VGammaDeexcitation* deexcitation = new G4DiscreteGammaDeexcitation;
    G4PhotonEvaporation* deexcitation = new G4PhotonEvaporation;
    deexcitation->SetVerboseLevel(GetVerboseLevel());
    deexcitation->Initialize(nucleus);
  //
  //
  // Get the gammas by deexciting the nucleus.
  //
    G4FragmentVector* gammas = deexcitation->BreakItUp(nucleus);
    G4int nGammas=0;
    //    if (gammas!=0) nGammas=gammas->entries();
    // in the case of BreakItUp(nucleus), the returned G4FragmentVector contains the residual nuclide
    // as its last entry.
    if (gammas > 0) nGammas=gammas->entries()-1;
  //
  //
  // Go through each gamma and add it to the decay product.  The angular distribution
  // of the gammas is isotropic, and the residual nucleus is assumed not to suffer
  // any recoil as a result of this de-excitation.
  //
    for (G4int ig=0; ig<nGammas; ig++)
    {
      G4double costheta = 2.0*G4UniformRand() - 1.0;
      G4double sintheta = sqrt((1.0 - costheta) * (1.0+costheta));
      G4double phi      = twopi * G4UniformRand();
      G4ParticleMomentum gDirection
        (sintheta*cos(phi),sintheta*sin(phi),costheta);
      G4double gEnergy = gammas->at(ig)->GetMomentum().e();
      G4DynamicParticle *theGammaRay = new
        G4DynamicParticle (G4Gamma::GammaDefinition(), gDirection, gEnergy);
      products->PushProducts (theGammaRay);
    }
  //
  //
  // Delete/reset variables associated with the gammas.
  //
    //    if (nGammas != 0) gammas->clearAndDestroy();
    gammas->clearAndDestroy();
    delete gammas;
    delete deexcitation;
  //
  //
  // Finally, add the residual nucleus to the product vector.  Note that for the
  // moment the daughter is assumed to end-up in the GROUND STATE, since
  // G4PhotoEvaporation cannot be made to stop when it reaches a long-lived
  // isomer.
  //
    G4double finalDaughterExcitation = 0.0*MeV;
    G4IonTable *theIonTable =  (G4IonTable*)(G4ParticleTable::GetParticleTable()->GetIonTable());
    dynamicDaughter = new G4DynamicParticle
      (theIonTable->GetIon(daughterZ,daughterA,finalDaughterExcitation),
      daughterMomentum1);
    products->PushProducts (dynamicDaughter); 
  }
  return products;
}
////////////////////////////////////////////////////////////////////////////////
//


G4DecayProducts *G4NuclearDecayChannel::BetaDecayIt()
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

  
  //calculate daughter momentum
  //  Generate two
  G4double rd1, rd2, rd;
  G4double daughtermomentum[3];
  G4double daughterenergy[3];
  G4double momentummax=0.0, momentumsum = 0.0;
  G4double Q = pmass - sumofdaughtermass;
  G4double fermif;
  G4BetaFermiFunction* aBetaFermiFunction;
  if (decayMode == 1) {
    // beta-decay 
    aBetaFermiFunction = new G4BetaFermiFunction (daughterA, daughterZ);
  } else {
    // beta+decay
    aBetaFermiFunction = new G4BetaFermiFunction (daughterA, -daughterZ);
  }
  if (GetVerboseLevel()>1) {
    G4cout<< " Q =  " <<Q<<G4endl;
    G4cout<< " daughterA =  " <<daughterA<<G4endl;
    G4cout<< " daughterZ =  " <<daughterZ<<G4endl;
    G4cout<< " decayMode = " <<decayMode << G4endl;
    G4cout<< " FermiFN =  " <<FermiFN<<G4endl;
  }
  do
    {
     rd1 = G4UniformRand();
     rd2 = G4UniformRand();

     momentummax = 0.0;
     momentumsum = 0.0;

     // daughter 0

     //     energy = rd2*(pmass - sumofdaughtermass);
     daughtermomentum[0] = sqrt(rd2) * sqrt((Q + 2.0*daughtermass[0])*Q);
     daughterenergy[0] = sqrt(daughtermomentum[0]*daughtermomentum[0] +
                         daughtermass[0] * daughtermass[0]) - daughtermass[0];
     if ( daughtermomentum[0] >momentummax )momentummax =  daughtermomentum[0];
     momentumsum  +=  daughtermomentum[0];

     // daughter 2
     //     energy = (1.-rd1)*(pmass - sumofdaughtermass);
     daughtermomentum[2] = sqrt(rd1)*sqrt((Q + 2.0*daughtermass[2])*Q);
     daughterenergy[2] = sqrt(daughtermomentum[2]*daughtermomentum[2] +
                         daughtermass[2] * daughtermass[2]) - daughtermass[2];
     if ( daughtermomentum[2] >momentummax )momentummax =  daughtermomentum[2];
     momentumsum  +=  daughtermomentum[2];

     // daughter 1

     daughterenergy[1] = Q - daughterenergy[0] - daughterenergy[2];
     if (daughterenergy[1] > 0.0) {
         daughtermomentum[1] = sqrt(daughterenergy[1]*daughterenergy[1] +
                               2.0*daughterenergy[1] * daughtermass[1]);
	 if ( daughtermomentum[1] >momentummax ) momentummax =
            daughtermomentum[1];
	 momentumsum +=  daughtermomentum[1];
     } else {
       momentummax = momentumsum = Q;
     }
     // beta particles is sampled with no coulomb effects applied above. Now
     // apply the Fermi function using rejection method.
     if (FermiOn) {
       daughterenergy[0] = daughterenergy[0]*MeV/0.511;
       fermif = aBetaFermiFunction->GetFF(daughterenergy[0])/FermiFN;
       // fermif: normalised Fermi factor
       if (G4UniformRand() > fermif) momentummax = momentumsum =  Q;
       // rejection method
     }
    } while (momentummax >  momentumsum - momentummax );

  delete aBetaFermiFunction;
  
  // output message
  if (GetVerboseLevel()>1) {
    G4cout <<"     daughter 0:" <<daughtermomentum[0]/GeV <<"[GeV/c]" <<G4endl;
    G4cout <<"     daughter 1:" <<daughtermomentum[1]/GeV <<"[GeV/c]" <<G4endl;
    G4cout <<"     daughter 2:" <<daughtermomentum[2]/GeV <<"[GeV/c]" <<G4endl;
    G4cout <<"   momentum sum:" <<momentumsum/GeV <<"[GeV/c]" <<G4endl;
  }

  //create daughter G4DynamicParticle
  G4double costheta, sintheta, phi, sinphi, cosphi;
  G4double costhetan, sinthetan, phin, sinphin, cosphin;
  costheta = 2.*G4UniformRand()-1.0;
  sintheta = sqrt((1.0-costheta)*(1.0+costheta));
  phi  = 2.0*M_PI*G4UniformRand()*rad;
  sinphi = sin(phi);
  cosphi = cos(phi);
  G4ParticleMomentum direction0(sintheta*cosphi,sintheta*sinphi,costheta);
  G4DynamicParticle * daughterparticle
         = new G4DynamicParticle( daughters[0], direction0*daughtermomentum[0]);
  products->PushProducts(daughterparticle);

  costhetan = (daughtermomentum[1]*daughtermomentum[1]-
    daughtermomentum[2]*daughtermomentum[2]-
    daughtermomentum[0]*daughtermomentum[0])/
    (2.0*daughtermomentum[2]*daughtermomentum[0]);
  sinthetan = sqrt((1.0-costhetan)*(1.0+costhetan));
  phin  = 2.0*M_PI*G4UniformRand()*rad;
  sinphin = sin(phin);
  cosphin = cos(phin);
  G4ParticleMomentum direction2;
  direction2.setX( sinthetan*cosphin*costheta*cosphi - 
    sinthetan*sinphin*sinphi + costhetan*sintheta*cosphi);
  direction2.setY( sinthetan*cosphin*costheta*sinphi +
    sinthetan*sinphin*cosphi + costhetan*sintheta*sinphi);
  direction2.setZ( -sinthetan*cosphin*sintheta +
    costhetan*costheta);
  daughterparticle = new G4DynamicParticle
    ( daughters[2], direction2*(daughtermomentum[2]/direction2.mag()));
  products->PushProducts(daughterparticle);

  daughterparticle =
       new G4DynamicParticle (daughters[1],
              (direction0*daughtermomentum[0] +
               direction2*(daughtermomentum[2]/direction2.mag()))*(-1.0));
  products->PushProducts(daughterparticle);

  //   delete daughterparticle;
  
  if (GetVerboseLevel()>1) {
     G4cout << "G4NuclearDecayChannel::BetaDecayIt ";
     G4cout << "  create decay products in rest frame " <<G4endl;
     products->DumpInfo();
  }
  return products;}
  








