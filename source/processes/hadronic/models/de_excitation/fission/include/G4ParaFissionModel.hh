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
#ifndef G4ParaFissionModel_h
#define G4ParaFissionModel_h 1

#include "G4CompetitiveFission.hh"
#include "G4ExcitationHandler.hh"
#include "G4HadronicInteraction.hh"
#include "G4NucleiProperties.hh"
//#include "G4ParticleTable.hh"

// Class Description
// Final state production model for (based on evaluated data
// libraries) description of neutron induced fission below 60 MeV; 
// In case you need the fission fragments, use this model.
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.


class G4ParaFissionModel : public G4HadronicInteraction
{
public:

  G4ParaFissionModel() 
  {
    SetMinEnergy( 0.0 );
    SetMaxEnergy( 60.*MeV );
  }

  virtual ~G4ParaFissionModel() {};
  
  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
					 G4Nucleus& theNucleus)
  {
    theParticleChange.Clear();
    theParticleChange.SetStatusChange( stopAndKill );
    theParticleChange.SetEnergyChange( 0.0 );
    
    // prepare the fragment

    G4int A = theNucleus.GetA_asInt();
    G4int Z = theNucleus.GetZ_asInt();
    G4double nucMass = G4NucleiProperties::GetNuclearMass(A, Z);
     
    G4int numberOfEx = aTrack.GetDefinition()->GetBaryonNumber();
    G4int numberOfCh = G4int(aTrack.GetDefinition()->GetPDGCharge() + 0.5);
    G4int numberOfHoles = 0;

    A += numberOfEx;
    Z += numberOfCh;
     
    G4LorentzVector v = aTrack.Get4Momentum() + G4LorentzVector(0.0,0.0,0.0,nucMass);
    G4Fragment anInitialState(A,Z,v);
    anInitialState.SetNumberOfExcitedParticle(numberOfEx,numberOfCh); 
    anInitialState.SetNumberOfHoles(0,0);

    // do the fission
    G4FragmentVector * theFissionResult = theFission.BreakUp(anInitialState);
    
    // deexcite the fission fragments and fill result

    G4int ll = theFissionResult->size();
    for(G4int i=0; i<ll; i++)
    {
      G4ReactionProductVector* theExcitationResult = 0; 
      G4Fragment* aFragment = (*theFissionResult)[i];
      if(aFragment->GetExcitationEnergy() > keV)
      {
	theExcitationResult = theHandler.BreakItUp(*aFragment);

	// add secondaries
	for(G4int j = 0; j < G4int(theExcitationResult->size()); j++)
	{
	  G4ReactionProduct* rp0 = (*theExcitationResult)[j];
	  G4DynamicParticle* p0 = 
	    new G4DynamicParticle(rp0->GetDefinition(),rp0->GetMomentum());
	  theParticleChange.AddSecondary(p0);
	  delete rp0;
	}
	delete theExcitationResult;
      }
      else
      {
	// add secondary
	G4DynamicParticle* p0 = 
	  new G4DynamicParticle(aFragment->GetParticleDefinition(),
				aFragment->GetMomentum());
	theParticleChange.AddSecondary(p0);
      }
      delete aFragment;
    }

    delete theFissionResult;
    
    return &theParticleChange;
  }
private:

  G4CompetitiveFission theFission;
  G4ExcitationHandler theHandler;
  
  G4HadFinalState theParticleChange;
};
#endif
