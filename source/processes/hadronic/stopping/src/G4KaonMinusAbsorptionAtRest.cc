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
//      Author:        Christian V"olcker (Christian.Volcker@cern.ch),
// 
//      Creation date: November 1997
//
//      Testfile:     ../G4KaonMinusAbsorptionAtRestTest.cc
//       
//      Modifications: 
//      Maria Grazia Pia  September 1998
//                  Various bug fixes, eliminated several memory leaks
//
// -------------------------------------------------------------------


#include "G4KaonMinusAbsorptionAtRest.hh"

#include "G4StopDeexcitation.hh"
#include "G4StopTheoDeexcitation.hh"
#include "G4StopDeexcitationAlgorithm.hh"
#include "G4ReactionKinematics.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicDeprecate.hh"

G4KaonMinusAbsorptionAtRest::G4KaonMinusAbsorptionAtRest(const G4String& processName,
                                      G4ProcessType   aType ) :
  G4VRestProcess (processName, aType)
{
  G4HadronicDeprecate("G4KaonMinusAbsorptionAtRest");
  if (verboseLevel>0) {
    G4cout << GetProcessName() << " is created "<< G4endl;
  }
  SetProcessSubType(fHadronAtRest);

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

  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}


G4KaonMinusAbsorptionAtRest::~G4KaonMinusAbsorptionAtRest()
{ 
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
}

void G4KaonMinusAbsorptionAtRest::PreparePhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

void G4KaonMinusAbsorptionAtRest::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

G4VParticleChange* G4KaonMinusAbsorptionAtRest::AtRestDoIt
(const G4Track& track, const G4Step& )
{
  stoppedHadron = track.GetDynamicParticle();
  
  // Check applicability

  if (!IsApplicable(*(stoppedHadron->GetDefinition()))) 
    {
      G4cerr  <<"G4KaonMinusAbsorptionAtRest:ERROR, particle must be a Kaon!" <<G4endl;
      return 0;
    }
  
  G4Material* material;
  material = track.GetMaterial();
  nucleus = 0;
  do
    {
      // Select the nucleus, get nucleon
      nucleus = new G4Nucleus(material);
      if (nucleus->GetA_asInt() < 1.5)
        {
          delete nucleus;
          nucleus = 0;
        }
    }  while(nucleus == 0);
    
  G4double Z = nucleus->GetZ_asInt();
  G4double A = nucleus->GetA_asInt();

  // Do the interaction with the nucleon
  G4DynamicParticleVector* absorptionProducts = KaonNucleonReaction();

  //A.R. 26-Jul-2012 Coverity fix
  if ( ! absorptionProducts ) {
    G4Exception("G4KaonMinusAbsorptionAtRest::AtRestDoIt()", "HAD_STOP_0001",
                FatalException, "NULL absorptionProducts");
    return 0;
  }
  
  // Secondary interactions
  
  G4DynamicParticle* thePion;
  unsigned int i;
  for(i = 0; i < absorptionProducts->size(); i++)
    {
      thePion = (*absorptionProducts)[i];
      if (thePion->GetDefinition() == G4PionMinus::PionMinus()
	  || thePion->GetDefinition() == G4PionPlus::PionPlus()
	  || thePion->GetDefinition() == G4PionZero::PionZero()) 
	{
	  if (AbsorbPionByNucleus(thePion))
	    {
	      absorptionProducts->erase(absorptionProducts->begin()+i);
	      i--;
              delete thePion;
	      if (verboseLevel > 1) 
		G4cout << "G4KaonMinusAbsorption::AtRestDoIt: Pion absorbed in Nucleus" 
		       << G4endl;
	    }                 
	}
    }
  
  G4DynamicParticle* theSigma;
  G4DynamicParticle* theLambda;
  for (i = 0; i < absorptionProducts->size(); i++)
    {
      theSigma = (*absorptionProducts)[i];
      if (theSigma->GetDefinition() == G4SigmaMinus::SigmaMinus()
	  || theSigma->GetDefinition() == G4SigmaPlus::SigmaPlus()
	  || theSigma->GetDefinition() == G4SigmaZero::SigmaZero()) 
	{
	  theLambda = SigmaLambdaConversion(theSigma);
	  if (theLambda  != 0){
	    absorptionProducts->erase(absorptionProducts->begin()+i);
	    i--;
            delete theSigma;
	    absorptionProducts->push_back(theLambda);

	    if (verboseLevel > 1) 
	      G4cout << "G4KaonMinusAbsorption::AtRestDoIt: SigmaLambdaConversion Done" 
		     << G4endl;
	  }                 
	}
    }
  
  // Nucleus deexcitation
  
  G4double productEnergy = 0.;
  G4ThreeVector pProducts(0.,0.,0.);

  unsigned int nAbsorptionProducts = 0;
  if (absorptionProducts != 0) nAbsorptionProducts = absorptionProducts->size();
  
  for ( i = 0; i<nAbsorptionProducts; i++)
    {
      pProducts += (*absorptionProducts)[i]->GetMomentum();
      productEnergy += (*absorptionProducts)[i]->GetKineticEnergy();
    }

  G4double newZ = nucleus->GetZ_asInt();
  G4double newA = nucleus->GetA_asInt();

  G4double bDiff = G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(A),static_cast<G4int>(Z)) - 
    G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(newA), static_cast<G4int>(newZ));
  
  G4StopDeexcitationAlgorithm* nucleusAlgorithm = new G4StopTheoDeexcitation();
  G4StopDeexcitation stopDeexcitation(nucleusAlgorithm);

  nucleus->AddExcitationEnergy(bDiff);
   
  // returns excitation energy for the moment ..
  G4double energyDeposit = nucleus->GetEnergyDeposit(); 
  if (verboseLevel>0)  
    {
      G4cout << " -- KaonAtRest -- excitation = " 
	     << energyDeposit 
	     << ", pNucleus = "
	     << pProducts 
	     << ", A: "
	     << A 
	     << ", "
	     << newA 
	     << ", Z: "
	     << Z
	     << ", "
	     << newZ
	     << G4endl; 
    }

  if (energyDeposit < 0.)
    G4Exception("G4KaonMinusAbsorptionAtRest::AtRestDoIt()", "HAD_STOP_0002",
                FatalException, "Excitation energy < 0");
  delete nucleus;    

  G4ReactionProductVector* fragmentationProducts = stopDeexcitation.DoBreakUp(newA,newZ,energyDeposit,pProducts);
  
  unsigned int nFragmentationProducts = 0;
  if (fragmentationProducts != 0) nFragmentationProducts = fragmentationProducts->size();
  
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
       new G4DynamicParticle((*fragmentationProducts)[i]->GetDefinition(),
                             (*fragmentationProducts)[i]->GetTotalEnergy(),
                             (*fragmentationProducts)[i]->GetMomentum());
    G4double newTime = aParticleChange.GetGlobalTime((*fragmentationProducts)[i]->GetFormationTime());
    aParticleChange.AddSecondary(aNew, newTime);
    delete (*fragmentationProducts)[i];
  }
  if (fragmentationProducts != 0) delete fragmentationProducts;
  
  // finally ...
  aParticleChange.ProposeTrackStatus(fStopAndKill); // Kill the incident Kaon
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
  
  G4double myZ = nucleus->GetZ_asInt();   // number of protons
  G4double myN = nucleus->GetA_asInt();   // number of nucleons (not neutrons!!)
  
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
  
  G4double iniZ = nucleus->GetZ_asInt();
  G4double iniA = nucleus->GetA_asInt();   
  
  G4DynamicParticle aNucleon = GetAbsorbingNucleon();
  
  // DHW 15 may 2011: unused: G4double nucleonMass;
  
  if (aNucleon.GetDefinition() == G4Proton::Proton()) 
    {
      // DHW 15 May 2011: unused: nucleonMass = proton_mass_c2+electron_mass_c2;
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
      // DHW 15 May 2011: unused: nucleonMass = neutron_mass_c2;
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
	    << G4endl;
	}

      //A.R. 26-Jul-2012 Coverity fix
      if ( products ) delete products;

      return 0;
    }  

  G4double newZ = nucleus->GetZ_asInt();
  G4double newA = nucleus->GetA_asInt();   
  
  // Modify the Kaon mass to take nuclear binding energy into account  
  // .. using mas formula ..
  // .. using mass table ..
  // equivalent to -'initialBindingEnergy+nucleus.GetBindingEnergy' !

  G4double nucleonBindingEnergy = 
    -G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(iniA), static_cast<G4int>(iniZ) )
    +G4NucleiProperties::GetBindingEnergy(static_cast<G4int>(newA), static_cast<G4int>(newZ) );
  
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
  
  products->push_back(producedBaryon);
  products->push_back(producedMeson);
  
  if (verboseLevel > 1) 
    {
      G4cout 
	<< "G4KaonMinusAbsorption::KaonNucleonReaction: Number of primaries = " 
	<< products->size()
	<< ": " <<producedMesonDef->GetParticleName() 
	<< ", " <<producedBaryonDef->GetParticleName() << G4endl;
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
  
  G4double A = nucleus->GetA_asInt();
  G4double Z = nucleus->GetZ_asInt();
  
  G4double newZ = Z;
  // DHW 15 May 2011: unused: G4double nucleonMassDifference = 0;
  
  G4ParticleDefinition* inNucleonDef=NULL;
  G4ParticleDefinition* outNucleonDef=NULL;

  // Decide which sigma
  switch((int) aSigma->GetDefinition()->GetPDGCharge()) {

  case 1: 
    sigmaLambdaConversionRate = sigmaPlusLambdaConversionRate;
    inNucleonDef   = G4Neutron::Neutron();
    outNucleonDef  = G4Proton::Proton();
    newZ = Z+1;
    // DHW 15 May 2011: unused: nucleonMassDifference =   neutron_mass_c2 - proton_mass_c2-electron_mass_c2;
    break;

  case -1: 
    sigmaLambdaConversionRate = sigmaMinusLambdaConversionRate;
    inNucleonDef   = G4Proton::Proton();
    outNucleonDef  = G4Neutron::Neutron();
    newZ = Z-1;
    // DHW 15 May 2011: unused: nucleonMassDifference =  proton_mass_c2+electron_mass_c2 - neutron_mass_c2;
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
    // Add dummy particles to avoid possibility of passing NULL pointers
    inNucleonDef   = G4Proton::Proton();
    outNucleonDef  = G4Proton::Proton();
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
  // .. using mass table ..
  // equivalent to -'initialBindingEnergy+nucleus.GetBindingEnergy' !

  // Add energy and momentum to nucleus, change Z,A 
  nucleus->AddExcitationEnergy(outNucleon.GetKineticEnergy());
  nucleus->AddMomentum(outNucleon.GetMomentum());
  nucleus->SetParameters(A,newZ);
  
  // The calling routine is responsible to delete the sigma!!
  return outLambda;
}



