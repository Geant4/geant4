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
// $Id$
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

G4TheoFSGenerator::G4TheoFSGenerator(const G4String& name)
    : G4HadronicInteraction(name)
    , theTransport(0), theHighEnergyGenerator(0)
    , theQuasielastic(0)
 {
 theParticleChange = new G4HadFinalState;
}

G4TheoFSGenerator::~G4TheoFSGenerator()
{
  delete theParticleChange;
}

void G4TheoFSGenerator::ModelDescription(std::ostream& outFile) const
{
  outFile << GetModelName() <<" consists of a " << theHighEnergyGenerator->GetModelName()
		  << " string model and a stage to de-excite the excited nuclear fragment."
		  << "\n<p>"
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
  
  // check if models have been registered, and use default, in case this is not true @@
  
  const G4DynamicParticle aPart(thePrimary.GetDefinition(),thePrimary.Get4Momentum().vect());

  if ( theQuasielastic ) {
  
     if ( theQuasielastic->GetFraction(theNucleus, aPart) > G4UniformRand() )
     {
       //G4cout<<"___G4TheoFSGenerator: before Scatter (1) QE=" << theQuasielastic<<G4endl;
       G4KineticTrackVector *result= theQuasielastic->Scatter(theNucleus, aPart);
       //G4cout << "^^G4TheoFSGenerator: after Scatter (1) " << G4endl;
       if (result)
       {
	    for(unsigned int  i=0; i<result->size(); i++)
	    {
	      G4DynamicParticle * aNew = 
		 new G4DynamicParticle(result->operator[](i)->GetDefinition(),
                        	       result->operator[](i)->Get4Momentum().e(),
                        	       result->operator[](i)->Get4Momentum().vect());
	      theParticleChange->AddSecondary(aNew);
	      delete result->operator[](i);
	    }
	    delete result;
	   
       } else 
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

//#define DEBUG_initial_result
  #ifdef DEBUG_initial_result
  	  G4double E_out(0);
  	  G4IonTable * ionTable=G4ParticleTable::GetParticleTable()->GetIonTable();
  	  std::vector<G4KineticTrack *>::iterator ir_iter;
  	  for(ir_iter=theInitialResult->begin(); ir_iter!=theInitialResult->end(); ir_iter++)
  	  {
  		  //G4cout << "TheoFS secondary, mom " << (*ir_iter)->GetDefinition()->GetParticleName() << " " << (*ir_iter)->Get4Momentum() << G4endl;
  		  E_out += (*ir_iter)->Get4Momentum().e();
  	  }
  	  G4double init_mass= ionTable->GetIonMass(theNucleus.GetZ_asInt(),theNucleus.GetA_asInt());
          G4double init_E=aPart.Get4Momentum().e();
  	  // residual nucleus

  	  const std::vector<G4Nucleon> & thy = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();

  	  G4int resZ(0),resA(0);
	  G4double delta_m(0);
  	  for(size_t them=0; them<thy.size(); them++)
  	  {
   	     if(thy[them].AreYouHit()) {
  	       ++resA;
  	       if ( thy[them].GetDefinition() == G4Proton::Proton() ) {
	          ++resZ;
		  delta_m +=G4Proton::Proton()->GetPDGMass();
	       } else {
	          delta_m +=G4Neutron::Neutron()->GetPDGMass();
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

  G4ReactionProductVector * theTransportResult = NULL;

// Uzhi Nov. 2012
  G4V3DNucleus* theProjectileNucleus = theHighEnergyGenerator->GetProjectileNucleus(); 
if(theProjectileNucleus == 0)                                       // Uzhi Nov. 2012
{                                                                   // Uzhi Nov. 2012

  G4int hitCount = 0;
  const std::vector<G4Nucleon>& they = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();
  for(size_t them=0; them<they.size(); them++)
  {
    if(they[them].AreYouHit()) hitCount ++;
  }
  if(hitCount != theHighEnergyGenerator->GetWoundedNucleus()->GetMassNumber() )
  {
    theTransport->SetPrimaryProjectile(thePrimary);	// For Bertini Cascade
    theTransportResult = 
               theTransport->Propagate(theInitialResult, theHighEnergyGenerator->GetWoundedNucleus());
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from transport propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");
    } 
  }
  else
  {
    theTransportResult = theDecay.Propagate(theInitialResult, theHighEnergyGenerator->GetWoundedNucleus());
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from decay propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from decay propagate");
    }   
  }

} else                                                              // Uzhi Nov. 2012
{                                                                   // Uzhi Nov. 2012
    theTransport->SetPrimaryProjectile(thePrimary);
    theTransportResult = 
    theTransport->PropagateNuclNucl(theInitialResult, 
                            theHighEnergyGenerator->GetWoundedNucleus(),
                            theHighEnergyGenerator->GetProjectileNucleus());
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from transport propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");
    } 
}                                                                   // Uzhi Nov. 2012

  // Fill particle change
  for(auto i=theTransportResult->begin(); i!=theTransportResult->end(); i++)
  {
    G4DynamicParticle * aNewDP =
       new G4DynamicParticle((*i)->GetDefinition(),
                             (*i)->GetTotalEnergy(),
                             (*i)->GetMomentum());
	G4HadSecondary aNew = G4HadSecondary(aNewDP);
    G4double time=(*i)->GetFormationTime();
    if(time < 0.0) { time = 0.0; }
    aNew.SetTime(timePrimary + time);
    aNew.SetCreatorModelType((*i)->GetCreatorModel());
    theParticleChange->AddSecondary(aNew);
    delete (*i);
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
