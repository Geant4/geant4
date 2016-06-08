// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//
// -------------------------------------------------------------------
//      GEANT 4 class file --- Copyright CERN 1997
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//
//      File name:     G4KaonMinusAbsorptionAtRest.cc 
//
//      Author:        Christian V"olcker (Christian.Volcker@cern.ch),
// 
//      Creation date: November 1997
//
//      Testfile:     ../G4KaonMinusAbsorptionAtRestTest.cc
//       
//      Modifications: 
//      Maria Grazia Pia  September 1998
//                        Various bug fixes, eliminated several memory leaks
//
// -------------------------------------------------------------------


#include "G4KaonMinusAbsorptionAtRest.hh"

#include "G4StopDeexcitation.hh"
#include "G4StopTheoDeexcitation.hh"
#include "G4StopDeexcitationAlgorithm.hh"
#include "G4ReactionKinematics.hh"

G4KaonMinusAbsorptionAtRest::G4KaonMinusAbsorptionAtRest(const G4String& processName)
  : G4VRestProcess (processName)
{
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created "<< endl;
  }

  // see Cohn et al, PLB27(1968) 527;
  //     Davis et al, PLB1(1967) 434; 
  
  pionAbsorptionRate = 0.07;
  
  // see  VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  // see  VanderVelde-Wilquet et al, Nuov.Cim.38A(1977)178;
  // see  VanderVelde-Wilquet et al, Nucl.Phys.A241(1975)511;
  // primary production rates ( for absorption on Carbon)
  // .. other elements are extrapolated by the halo factor.
  
  rateLambdaZeroPiZero = 0.052;
  rateSigmaMinusPiPlus = 0.199;
  rateSigmaPlusPiMinus = 0.446;
  rateSigmaZeroPiZero  = 0.303;
  rateLambdaZeroPiMinus = 0.568;
  rateSigmaZeroPiMinus  = 0.216;
  rateSigmaMinusPiZero  = 0.216;
  
  // for sigma- p -> lambda n
  //     sigma+ n -> lambda p
  //     sigma- n -> lambda 
  // all values compatible with 0.55 same literature as above.
  
  sigmaPlusLambdaConversionRate = 0.55; 
  sigmaMinusLambdaConversionRate = 0.55;
  sigmaZeroLambdaConversionRate = 0.55;
}


G4KaonMinusAbsorptionAtRest::~G4KaonMinusAbsorptionAtRest()
{ }


G4VParticleChange* G4KaonMinusAbsorptionAtRest::AtRestDoIt
(const G4Track& track, const G4Step& Step)
{
  stoppedHadron = track.GetDynamicParticle();
  
  // Check applicability

  if (!IsApplicable(*(stoppedHadron->GetDefinition()))) 
    {
      G4cerr  <<"G4KaonMinusAbsorptionAtRest:ERROR, particle must be a Kaon!" <<endl;
      return 0;
    }
  
  G4Material* material;
  material = track.GetMaterial();
  nucleus = 0;
  do
    {
      // Select the nucleus, get nucleon
      nucleus = new G4Nucleus(material);
      if (nucleus->GetN() < 1.5)
        {
          delete nucleus;
          nucleus = 0;
        }
    }  while(nucleus == 0);
    
  G4double Z = nucleus->GetZ();
  G4double A = nucleus->GetN();

  // Do the interaction with the nucleon
  G4DynamicParticleVector* absorptionProducts = KaonNucleonReaction();
  
  // Secondary interactions
  
  G4DynamicParticle* thePion;
  G4int i;
  for(i = 0; i < absorptionProducts->length(); i++)
    {
      thePion = (*absorptionProducts)[i];
      if (thePion->GetDefinition() == G4PionMinus::PionMinus()
	  || thePion->GetDefinition() == G4PionPlus::PionPlus()
	  || thePion->GetDefinition() == G4PionZero::PionZero()) 
	{
	  if (AbsorbPionByNucleus(thePion))
	    {
	      absorptionProducts->remove(thePion);
              delete thePion;
	      if (verboseLevel > 1) 
		G4cout << "G4KaonMinusAbsorption::AtRestDoIt: Pion absorbed in Nucleus" 
		       << endl;
	    }                 
	}
    }
  
  G4DynamicParticle* theSigma;
  G4DynamicParticle* theLambda;
  for (i = 0; i < absorptionProducts->length(); i++)
    {
      theSigma = (*absorptionProducts)[i];
      if (theSigma->GetDefinition() == G4SigmaMinus::SigmaMinus()
	  || theSigma->GetDefinition() == G4SigmaPlus::SigmaPlus()
	  || theSigma->GetDefinition() == G4SigmaZero::SigmaZero()) 
	{
	  theLambda = SigmaLambdaConversion(theSigma);
	  if (theLambda  != 0){
	    absorptionProducts->remove(theSigma);
            delete theSigma;
	    absorptionProducts->append(theLambda);

	    if (verboseLevel > 1) 
	      G4cout << "G4KaonMinusAbsorption::AtRestDoIt: SigmaLambdaConversion Done" 
		     << endl;
	  }                 
	}
    }
  
  // Nucleus deexcitation
  
  G4double productEnergy = 0.;
  G4ThreeVector pProducts(0.,0.,0.);
  G4int nN = 0;
  G4int nP = 0;

  G4int nAbsorptionProducts = 0;
  if (absorptionProducts != 0) nAbsorptionProducts = absorptionProducts->entries();
  
  for ( i = 0; i<nAbsorptionProducts; i++)
    {
      pProducts = pProducts + (*absorptionProducts)[i]->GetMomentum();
      productEnergy += (*absorptionProducts)[i]->GetKineticEnergy();
    }

  G4double newZ = nucleus->GetZ();
  G4double newA = nucleus->GetN();

  G4double bDiff = G4NucleiPropertiesTable::GetBindingEnergy(Z,A) - 
    G4NucleiPropertiesTable::GetBindingEnergy(newZ,newA);

  //  G4double mass = G4NucleiPropertiesTable::GetAtomicMass(newZ,newA);
  G4double pNucleus = pProducts.mag();
  
  G4StopDeexcitationAlgorithm* nucleusAlgorithm = new G4StopTheoDeexcitation();
  G4StopDeexcitation stopDeexcitation(nucleusAlgorithm);

  // G4double difference = G4KaonMinus::KaonMinus()->GetPDGMass() - productEnergy - bDiff;

  nucleus->AddExcitationEnergy(bDiff);
   
  // returns excitation energy for the moment ..
  G4double energyDeposit = nucleus->GetEnergyDeposit(); 
  if (verboseLevel>0)  
    {
      G4cout << " -- KaonAtRest -- excitation = " 
	     << energyDeposit 
	     << ", pNucleus = "
	     << pNucleus 
	     << ", A: "
	     << A 
	     << ", "
	     << newA 
	     << ", Z: "
	     << Z
	     << ", "
	     << newZ
	     << endl; 
    }

  if (energyDeposit < 0.) 
      G4Exception("G4KaonMinusAbsorptionAtRest::AtRestDoIt -- excitation energy < 0");

  delete nucleus;    

  G4ReactionProductVector* fragmentationProducts = stopDeexcitation.DoBreakUp(newA,newZ,energyDeposit,pNucleus);
  
  G4int nFragmentationProducts = 0;
  if (fragmentationProducts != 0) nFragmentationProducts = fragmentationProducts->entries();
  
  //Initialize ParticleChange
   aParticleChange.Initialize(track);
  aParticleChange.SetNumberOfSecondaries(G4int(nAbsorptionProducts+nFragmentationProducts) ); 
  
  // update List of alive particles. put energy deposit at the right place ...
  for (i = 0; i < nAbsorptionProducts; i++)
    {aParticleChange.AddSecondary((*absorptionProducts)[i]); }
  if (absorptionProducts != 0) delete absorptionProducts;
  
//  for (i = 0; i < nFragmentationProducts; i++)
//    { aParticleChange.AddSecondary(fragmentationProducts->at(i)); }
  for(i=0; i<nFragmentationProducts; i++)
  {
    G4DynamicParticle * aNew = 
       new G4DynamicParticle(fragmentationProducts->at(i)->GetDefinition(),
                             fragmentationProducts->at(i)->GetTotalEnergy(),
                             fragmentationProducts->at(i)->GetMomentum());
    G4double newTime = aParticleChange.GetGlobalTime(fragmentationProducts->at(i)->GetFormationTime());
    aParticleChange.AddSecondary(aNew, newTime);
    delete fragmentationProducts->at(i);
  }
  if (fragmentationProducts != 0) delete fragmentationProducts;
  
  // finally ...
  aParticleChange.SetStatusChange(fStopAndKill); // Kill the incident Kaon
  return &aParticleChange;
}


G4DynamicParticle G4KaonMinusAbsorptionAtRest::GetAbsorbingNucleon()
{
  G4DynamicParticle aNucleon;
  
  // Get nucleon definition, based on Z,N of current Nucleus
  aNucleon.SetDefinition(SelectAbsorbingNucleon());
  
  // Fermi momentum distribution in three dimensions
  G4ThreeVector pFermi = nucleus->GetFermiMomentum();
  aNucleon.SetMomentum(pFermi);
  
  return aNucleon;
}

G4ParticleDefinition* G4KaonMinusAbsorptionAtRest::SelectAbsorbingNucleon()
{
  // (Ch. Voelcker) extended from ReturnTargetParticle():
  // Choose a proton or a neutron as the absorbing particle,
  // taking weight into account!
  // Update nucleon's atomic numbers.
  
  G4ParticleDefinition* absorbingParticleDef;
  
  G4double ranflat = G4UniformRand();   
  
  G4double myZ = nucleus->GetZ();   // number of protons
  G4double myN = nucleus->GetN();   // number of nucleons (not neutrons!!)
  
  // See  VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  G4double carbonRatioNP = 0.18;  // (Rn/Rp)c, see page 544 
  
  G4double neutronProtonRatio = NeutronHaloFactor(myZ,myN)*carbonRatioNP*(myN-myZ)/myZ;
  G4double protonProbability = 1./(1.+neutronProtonRatio);
  
  if ( ranflat < protonProbability ) 
    {
      absorbingParticleDef = G4Proton::Proton();
      myZ-= 1.;
    } 
  else 
    { absorbingParticleDef = G4Neutron::Neutron(); }

  myN -= 1.;
  nucleus->SetParameters(myN,myZ);
  return absorbingParticleDef;
}


G4double G4KaonMinusAbsorptionAtRest::NeutronHaloFactor(G4double Z, G4double N)
{
  // this function should take care of the probability for absorption
  // on neutrons, depending on number of protons Z and number of neutrons N-Z
  // parametrisation from fit to 
  // VanderVelde-Wilquet et al, Nuov.Cim.39A(1978)538;
  // 
  
  if (Z == 1.) return 1.389;      // deuterium
  else if (Z == 2.) return 1.78;  // helium
  else if (Z == 10.) return 0.66; // neon
  else     
    return 0.6742+(N-Z)*0.06524;  
}    


G4DynamicParticleVector* G4KaonMinusAbsorptionAtRest::KaonNucleonReaction()
{
  G4DynamicParticleVector* products = new G4DynamicParticleVector();
  
  G4double ranflat = G4UniformRand();   
  G4double prob = 0;
  
  G4ParticleDefinition* producedBaryonDef;
  G4ParticleDefinition* producedMesonDef;
  
  G4double iniZ = nucleus->GetZ();
  G4double iniA = nucleus->GetN();   
  
  G4DynamicParticle aNucleon = GetAbsorbingNucleon();
  
  G4double nucleonMass;
  
  if (aNucleon.GetDefinition() == G4Proton::Proton()) 
    {
      nucleonMass = proton_mass_c2+electron_mass_c2;
      if ( (prob += rateLambdaZeroPiZero) > ranflat) 
	{                                                  //  lambda pi0
	  producedBaryonDef = G4Lambda::Lambda();
	  producedMesonDef  = G4PionZero::PionZero();
	} 
      else if ((prob += rateSigmaPlusPiMinus) > ranflat) 
	{                                                  //  sigma+ pi-
	  producedBaryonDef = G4SigmaPlus::SigmaPlus();
	  producedMesonDef  = G4PionMinus::PionMinus();
	} 
      else if ((prob += rateSigmaMinusPiPlus) > ranflat) 
	{                                                  //  sigma- pi+
	  producedBaryonDef = G4SigmaMinus::SigmaMinus();
	  producedMesonDef  = G4PionPlus::PionPlus();
	} 
      else 
	{                                                 //  sigma0 pi0
	  producedBaryonDef = G4SigmaZero::SigmaZero();
	  producedMesonDef  = G4PionZero::PionZero();
	}
    } 
  else if (aNucleon.GetDefinition() == G4Neutron::Neutron()) 
    {
      nucleonMass = neutron_mass_c2;
      if ((prob += rateLambdaZeroPiMinus) > ranflat) 
	{                                                 //  lambda pi-
	  producedBaryonDef = G4Lambda::Lambda();
	  producedMesonDef  = G4PionMinus::PionMinus();
	} 
      else if ((prob += rateSigmaZeroPiMinus) > ranflat) 
	{                                                //  sigma0 pi-
	  producedBaryonDef = G4SigmaZero::SigmaZero();
	  producedMesonDef = G4PionMinus::PionMinus();
	} 
      else 
	{                                               //  sigma- pi0
	  producedBaryonDef = G4SigmaMinus::SigmaMinus();
	  producedMesonDef  = G4PionZero::PionZero();
	}
    } 
  else 
    {
      if (verboseLevel>0)
	{
	  G4cout 
	    << "G4KaonMinusAbsorption::KaonNucleonReaction: "
	    << aNucleon.GetDefinition()->GetParticleName() 
	    << " is not a good nucleon - check G4Nucleus::ReturnTargetParticle()!"
	    << endl;
	}
      return 0;
    }  

  G4double newZ = nucleus->GetZ();
  G4double newA = nucleus->GetN();   
  
  // Modify the Kaon mass to take nuclear binding energy into account  
  // .. using mas formula ..
  //   G4double nucleonBindingEnergy =  nucleus->AtomicMass(iniA,iniZ)
  //                                  - nucleus->AtomicMass(newA,newZ)
  //                                  - nucleonMass;
  // .. using mass table ..
  //   G4double nucleonBindingEnergy =  
  //            G4NucleiPropertiesTable::GetAtomicMass(iniZ,iniA)
  //           -G4NucleiPropertiesTable::GetAtomicMass(newZ,newA)
  //           -nucleonMass;
  // equivalent to -'initialBindingEnergy+nucleus.GetBindingEnergy' !

  G4double nucleonBindingEnergy = 
    -G4NucleiPropertiesTable::GetBindingEnergy(iniZ,iniA)
    +G4NucleiPropertiesTable::GetBindingEnergy(newZ,newA);
  
  G4DynamicParticle modifiedHadron = (*stoppedHadron);
  modifiedHadron.SetMass(stoppedHadron->GetMass() + nucleonBindingEnergy);   
  
  // Setup outgoing dynamic particles 
  G4ThreeVector dummy(0.,0.,0.);
  G4DynamicParticle* producedBaryon = new G4DynamicParticle(producedBaryonDef,dummy); 
  G4DynamicParticle* producedMeson = new G4DynamicParticle(producedMesonDef,dummy); 
  
  // Produce the secondary particles in a twobody process:
  G4ReactionKinematics theReactionKinematics;
  theReactionKinematics.TwoBodyScattering( &modifiedHadron, &aNucleon,
					   producedBaryon, producedMeson);
  
  products->append(producedBaryon);
  products->append(producedMeson);
  
  if (verboseLevel > 1) 
    {
      G4cout 
	<< "G4KaonMinusAbsorption::KaonNucleonReaction: Number of primaries = " 
	<< products->entries()
	<< ": " <<producedMesonDef->GetParticleName() 
	<< ", " <<producedBaryonDef->GetParticleName() << endl;
    }
  
  return products;
}


G4bool G4KaonMinusAbsorptionAtRest::AbsorbPionByNucleus(G4DynamicParticle* aPion)
{
  // Needs some more investigation!

  G4double ranflat = G4UniformRand();   

  if (ranflat < pionAbsorptionRate){
    // Add pion energy to ExcitationEnergy and NucleusMomentum
    nucleus->AddExcitationEnergy(aPion->GetTotalEnergy());
    nucleus->AddMomentum(aPion->GetMomentum());
  }

  return (ranflat < pionAbsorptionRate);
}

G4DynamicParticle* G4KaonMinusAbsorptionAtRest::SigmaLambdaConversion(G4DynamicParticle* aSigma)
{
  G4double  ranflat = G4UniformRand();
  G4double  sigmaLambdaConversionRate;
  
  G4double A = nucleus->GetN();
  G4double Z = nucleus->GetZ();
  
  G4double newZ = Z;
  G4double nucleonMassDifference = 0;
  
  G4ParticleDefinition* inNucleonDef;
  G4ParticleDefinition* outNucleonDef;

  // Decide which sigma
  switch((int) aSigma->GetDefinition()->GetPDGCharge()) {

  case 1: 
    sigmaLambdaConversionRate = sigmaPlusLambdaConversionRate;
    inNucleonDef   = G4Neutron::Neutron();
    outNucleonDef  = G4Proton::Proton();
    newZ = Z+1;
    nucleonMassDifference =   neutron_mass_c2 - proton_mass_c2-electron_mass_c2;
    break;

  case -1: 
    sigmaLambdaConversionRate = sigmaMinusLambdaConversionRate;
    inNucleonDef   = G4Proton::Proton();
    outNucleonDef  = G4Neutron::Neutron();
    newZ = Z-1;
    nucleonMassDifference =  proton_mass_c2+electron_mass_c2 - neutron_mass_c2;
    break;

  case 0: 
    sigmaLambdaConversionRate = sigmaZeroLambdaConversionRate;
    // The 'outgoing' nucleon is just virtual, to keep the energy-momentum 
    // balance and will not appear in the ParticleChange. Therefore no need 
    // choose between neutron and proton here!
    inNucleonDef   = G4Neutron::Neutron();
    outNucleonDef  = G4Neutron::Neutron();
    break;

  default: 
    sigmaLambdaConversionRate = 0.;
  }
  
  if (ranflat >= sigmaLambdaConversionRate) return 0;
  
  G4ThreeVector dummy(0.,0.,0.);
  
  // Fermi momentum distribution in three dimensions
  G4ThreeVector momentum = nucleus->GetFermiMomentum();

  G4ParticleDefinition* lambdaDef  = G4Lambda::Lambda();
  
  G4DynamicParticle inNucleon(inNucleonDef,momentum); 
  G4DynamicParticle outNucleon(outNucleonDef,dummy); 
  G4DynamicParticle* outLambda = new G4DynamicParticle(lambdaDef,dummy); 
  
  G4ReactionKinematics theReactionKinematics;
 
  // Now do the twobody scattering
  theReactionKinematics.TwoBodyScattering(aSigma, &inNucleon,
					  &outNucleon, outLambda);

  // Binding energy of nucleus has changed. This will change the
  // ExcitationEnergy.
  // .. using mass formula ..
  //   G4double massDifference =   nucleus->AtomicMass(A,Z)
  //                           - nucleus->AtomicMass(A,newZ)
  //                           - nucleonMassDifference;
  // .. using mass table ..
  // G4double massDifference = 
  //            G4NucleiPropertiesTable::GetAtomicMass(Z,A)
  //           -G4NucleiPropertiesTable::GetAtomicMass(newZ,A)
  //           -nucleonMass;
  // equivalent to -'initialBindingEnergy+nucleus.GetBindingEnergy' !
  G4double massDifference =
    -G4NucleiPropertiesTable::GetBindingEnergy(Z,A)
    +G4NucleiPropertiesTable::GetBindingEnergy(newZ,A);
  
  
  // Add energy and momentum to nucleus, change Z,A 
  //  nucleus->AddExcitationEnergy(outNucleon->GetKineticEnergy()+massDifference);
  nucleus->AddExcitationEnergy(outNucleon.GetKineticEnergy());
  nucleus->AddMomentum(outNucleon.GetMomentum());
  nucleus->SetParameters(A,newZ);
  
  // The calling routine is responsible to delete the sigma!!
  return outLambda;
}



