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
// G4TheoFSGenerator
//
// 20110307  M. Kelsey -- Add call to new theTransport->SetPrimaryProjectile()
//		to provide access to full initial state (for Bertini)
// 20110805  M. Kelsey -- Follow change to G4V3DNucleus::GetNucleons()

#include "G4DynamicParticle.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4IonTable.hh"
#include "G4HadronicParameters.hh"
#include "G4CRCoalescence.hh"
#include "G4HadronicInteractionRegistry.hh"

G4TheoFSGenerator::G4TheoFSGenerator(const G4String& name)
    : G4HadronicInteraction(name)
    , theTransport(nullptr), theHighEnergyGenerator(nullptr)
    , theQuasielastic(nullptr)
    , theCosmicCoalescence(nullptr)
    , theStringModelID(-1)
{
  theParticleChange = new G4HadFinalState;
  theStringModelID = G4PhysicsModelCatalog::GetModelID( "model_" + name );
}

G4TheoFSGenerator::~G4TheoFSGenerator()
{
  delete theParticleChange;
}

void G4TheoFSGenerator::ModelDescription(std::ostream& outFile) const
{
  outFile << GetModelName() <<" consists of a " << theHighEnergyGenerator->GetModelName()
	  << " string model and a stage to de-excite the excited nuclear fragment.\n<p>"
	  << "The string model simulates the interaction of\n"
          << "an incident hadron with a nucleus, forming \n"
          << "excited strings, decays these strings into hadrons,\n"
          << "and leaves an excited nucleus. \n"
          << "<p>The string model:\n";
  theHighEnergyGenerator->ModelDescription(outFile);
  outFile <<"\n<p>";
  theTransport->PropagateModelDescription(outFile);
}

G4HadFinalState * G4TheoFSGenerator::ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus &theNucleus)
{
  // init particle change
  theParticleChange->Clear();
  theParticleChange->SetStatusChange(stopAndKill);
  G4double timePrimary=thePrimary.GetGlobalTime();

  // Temporarily dummy treatment of heavy (charm and bottom) hadron projectiles at low energies.
  // Cascade models are currently not applicable for heavy hadrons and string models cannot
  // handle them properly at low energies - let's say safely below ~100 MeV.
  // In these cases, we return as final state the initial state unchanged.
  // For most applications, this is a safe simplification, giving that the nearly all
  // slowly moving charm and bottom hadrons decay before any hadronic interaction can occur.
  // Note that we prefer not to use G4HadronicParameters::GetMinEnergyTransitionFTF_Cascade()
  // (typically ~3 GeV) because FTFP works reasonably well below such a value.
  const G4double energyThresholdForCharmAndBottomHadrons = 100.0*CLHEP::MeV;
  if ( thePrimary.GetKineticEnergy() < energyThresholdForCharmAndBottomHadrons  &&
       ( thePrimary.GetDefinition()->GetQuarkContent( 4 )     != 0  ||    // Has charm       constituent quark
	 thePrimary.GetDefinition()->GetAntiQuarkContent( 4 ) != 0  ||    // Has anti-charm  constituent anti-quark
	 thePrimary.GetDefinition()->GetQuarkContent( 5 )     != 0  ||    // Has bottom      constituent quark
	 thePrimary.GetDefinition()->GetAntiQuarkContent( 5 ) != 0 ) ) {  // Has anti-bottom constituent anti-quark
    theParticleChange->SetStatusChange( isAlive );
    theParticleChange->SetEnergyChange( thePrimary.GetKineticEnergy() );
    theParticleChange->SetMomentumChange( thePrimary.Get4Momentum().vect().unit() );
    return theParticleChange;
  }

  // Temporarily dummy treatment of light hypernuclei projectiles at low energies.
  // Bertini and Binary intra-nuclear cascade models are currently not applicable
  // for hypernuclei and string models cannot handle them properly at low energies -
  // let's say safely below ~100 MeV.
  // In these cases, we return as final state the initial state unchanged.
  // For most applications, this is a safe simplification, giving that the nearly all
  // slowly moving hypernuclei decay before any hadronic interaction can occur.
  // Note that we prefer not to use G4HadronicParameters::GetMinEnergyTransitionFTF_Cascade()
  // (typically ~3 GeV) because FTFP works reasonably well below such a value.
  // Note that for anti-hypernuclei, FTF model works fine down to zero kinetic energy,
  // so there is no need of any extra, dummy treatment.
  const G4double energyThresholdForHyperNuclei = 100.0*CLHEP::MeV;
  if ( thePrimary.GetKineticEnergy() < energyThresholdForHyperNuclei  &&
       thePrimary.GetDefinition()->IsHypernucleus() ) {
    theParticleChange->SetStatusChange( isAlive );
    theParticleChange->SetEnergyChange( thePrimary.GetKineticEnergy() );
    theParticleChange->SetMomentumChange( thePrimary.Get4Momentum().vect().unit() );
    return theParticleChange;
  }

  // check if models have been registered, and use default, in case this is not true @@
  
  const G4DynamicParticle aPart(thePrimary.GetDefinition(),thePrimary.Get4Momentum().vect());

  if ( theQuasielastic ) 
  {
     if ( theQuasielastic->GetFraction(theNucleus, aPart) > G4UniformRand() )
     {
       //G4cout<<"___G4TheoFSGenerator: before Scatter (1) QE=" << theQuasielastic<<G4endl;
       G4KineticTrackVector *result= theQuasielastic->Scatter(theNucleus, aPart);
       //G4cout << "^^G4TheoFSGenerator: after Scatter (1) " << G4endl;
       if (result)
       {
	    for(auto & ptr : *result)
	    {
	      G4DynamicParticle * aNew = 
		 new G4DynamicParticle(ptr->GetDefinition(),
                        	       ptr->Get4Momentum().e(),
                        	       ptr->Get4Momentum().vect());
	      theParticleChange->AddSecondary(aNew, ptr->GetCreatorModelID());
	      delete ptr;
	    }
	    delete result;	   
       } 
       else 
       {
	    theParticleChange->SetStatusChange(isAlive);
	    theParticleChange->SetEnergyChange(thePrimary.GetKineticEnergy());
	    theParticleChange->SetMomentumChange(thePrimary.Get4Momentum().vect().unit());
       }
       return theParticleChange;
     } 
  }

  // get result from high energy model
  G4KineticTrackVector * theInitialResult =
               theHighEnergyGenerator->Scatter(theNucleus, aPart);

  // Assign the creator model ID
  for ( auto & ptr : *theInitialResult ) {
    ptr->SetCreatorModelID( theStringModelID );
  }

  //#define DEBUG_initial_result
  #ifdef DEBUG_initial_result
  	  G4double E_out(0);
  	  G4IonTable * ionTable=G4ParticleTable::GetParticleTable()->GetIonTable();
  	  for(auto & ptr : *theInitialResult)
  	  {
  	     //G4cout << "TheoFS secondary, mom " << ptr->GetDefinition()->GetParticleName() 
             //         << " " << ptr->Get4Momentum() << G4endl;
  	     E_out += ptr->Get4Momentum().e();
  	  }
  	  G4double init_mass= ionTable->GetIonMass(theNucleus.GetZ_asInt(),theNucleus.GetA_asInt());
          G4double init_E=aPart.Get4Momentum().e();
  	  // residual nucleus

  	  const std::vector<G4Nucleon> & thy = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();

  	  G4int resZ(0),resA(0);
	  G4double delta_m(0);
  	  for(auto & nuc : thy)
  	  {
   	     if(nuc.AreYouHit()) {
  	       ++resA;
  	       if ( nuc.GetDefinition() == G4Proton::Proton() ) {
	         ++resZ;
		 delta_m += CLHEP::proton_mass_c2;
	       } else {
		 delta_m += CLHEP::neutron_mass_c2;
	       }  
  	     }
	  }

  	  G4double final_mass(0);
	  if ( theNucleus.GetA_asInt() ) {
	   final_mass=ionTable->GetIonMass(theNucleus.GetZ_asInt()-resZ,theNucleus.GetA_asInt()- resA);
  	  }
	  G4double E_excit=init_mass + init_E - final_mass - E_out;
	  G4cout << " Corrected delta mass " << init_mass - final_mass - delta_m << G4endl;
  	  G4cout << "initial E, mass = " << init_E << ", " << init_mass << G4endl;
  	  G4cout << "  final E, mass = " << E_out <<", " << final_mass << "  excitation_E " << E_excit << G4endl;
  #endif

  G4ReactionProductVector * theTransportResult = nullptr;

  G4V3DNucleus* theProjectileNucleus = theHighEnergyGenerator->GetProjectileNucleus(); 
  if(theProjectileNucleus == nullptr)
  {
    G4int hitCount = 0;
    const std::vector<G4Nucleon>& they = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();
    for(auto & nuc : they)
    {
      if(nuc.AreYouHit()) ++hitCount;
    }
    if(hitCount != theHighEnergyGenerator->GetWoundedNucleus()->GetMassNumber() )
    {
      theTransport->SetPrimaryProjectile(thePrimary);
      theTransportResult = 
        theTransport->Propagate(theInitialResult, 
                                theHighEnergyGenerator->GetWoundedNucleus());
      if ( !theTransportResult ) {
        G4cout << "G4TheoFSGenerator: null ptr from transport propagate " << G4endl;
        throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");
      }
    }
    else
    {
      theTransportResult = theDecay.Propagate(theInitialResult, 
                                 theHighEnergyGenerator->GetWoundedNucleus());
      if ( theTransportResult == nullptr ) {
         G4cout << "G4TheoFSGenerator: null ptr from decay propagate " << G4endl;
         throw G4HadronicException(__FILE__, __LINE__, "Null ptr from decay propagate");
      }   
    }
  } 
  else 
  { 
    theTransport->SetPrimaryProjectile(thePrimary);
    theTransportResult = theTransport->PropagateNuclNucl(theInitialResult, 
                            theHighEnergyGenerator->GetWoundedNucleus(),
                            theProjectileNucleus);
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from transport propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");
    } 
  }

  // If enabled, apply the Cosmic Rays (CR) coalescence to the list of secondaries produced so far.
  // This algorithm can form deuterons and antideuterons by coalescence of, respectively,
  // proton-neutron and antiproton-antineutron pairs close in momentum space.
  // This can be useful in particular for Cosmic Ray applications.
  if ( G4HadronicParameters::Instance()->EnableCRCoalescence() ) {
    if(nullptr == theCosmicCoalescence) {
      theCosmicCoalescence = (G4CRCoalescence*)
        G4HadronicInteractionRegistry::Instance()->FindModel("G4CRCoalescence");
      if(nullptr == theCosmicCoalescence) { 
	theCosmicCoalescence = new G4CRCoalescence();
      }
    }
    theCosmicCoalescence->SetP0Coalescence( thePrimary, theHighEnergyGenerator->GetModelName() );
    theCosmicCoalescence->GenerateDeuterons( theTransportResult );
  }
    
  // Fill particle change
  for(auto & ptr : *theTransportResult)
  {
    G4DynamicParticle * aNewDP =
       new G4DynamicParticle(ptr->GetDefinition(),
                             ptr->GetTotalEnergy(),
                             ptr->GetMomentum());
    G4HadSecondary aNew = G4HadSecondary(aNewDP);
    G4double time = std::max(ptr->GetFormationTime(), 0.0);
    aNew.SetTime(timePrimary + time);
    aNew.SetCreatorModelID(ptr->GetCreatorModelID());
    aNew.SetParentResonanceDef(ptr->GetParentResonanceDef());
    aNew.SetParentResonanceID(ptr->GetParentResonanceID());    
    theParticleChange->AddSecondary(aNew);
    delete ptr;
  }
  
  // some garbage collection
  delete theTransportResult;
  
  // Done
  return theParticleChange;
}

std::pair<G4double, G4double> G4TheoFSGenerator::GetEnergyMomentumCheckLevels() const
{
  if ( theHighEnergyGenerator ) {
	 return theHighEnergyGenerator->GetEnergyMomentumCheckLevels();
  } else {
	 return std::pair<G4double, G4double>(DBL_MAX, DBL_MAX);
  }
}
