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
#include "globals.hh"
#include "G4GammaParticipants.hh"
#include "G4LorentzVector.hh"
#include "G4V3DNucleus.hh"
#include <utility>

// Class G4GammaParticipants 

G4VSplitableHadron* G4GammaParticipants::SelectInteractions(const G4ReactionProduct  &thePrimary) 
{
	// Check reaction threshold  - goes to CheckThreshold

	theProjectileSplitable = new G4QGSMSplitableHadron(thePrimary, TRUE);
        theProjectileSplitable->SetStatus(1);

	G4LorentzVector aPrimaryMomentum(thePrimary.GetMomentum(), thePrimary.GetTotalEnergy());
	G4LorentzVector aTargetNMomentum(0.,0.,0.,938.);
	if((!(aPrimaryMomentum.e()>-1)) && (!(aPrimaryMomentum.e()<1)) )
	{
		throw G4HadronicException(__FILE__, __LINE__,
				"G4GammaParticipants::SelectInteractions: primary nan energy.");
	}
	G4double S = (aPrimaryMomentum + aTargetNMomentum).mag2();
	G4double ThresholdMass = thePrimary.GetMass() + 938.;
	ModelMode = SOFT;

	if (sqr(ThresholdMass + ThresholdParameter) > S)
	{
		ModelMode = DIFFRACTIVE;
	}

	if (sqr(ThresholdMass + QGSMThreshold) > S)
	{
		ModelMode = DIFFRACTIVE;
	}

	std::for_each(theInteractions.begin(), theInteractions.end(), DeleteInteractionContent());
	theInteractions.clear();
	G4int totalCuts = 0;

	#ifdef debug_G4GammaParticipants
		G4double eK = thePrimary.GetKineticEnergy()/GeV;
		G4int nucleonCount = theNucleus->GetMassNumber();
	#endif

	G4int theCurrent = G4int(theNucleus->GetMassNumber()*G4UniformRand());
        G4int NucleonNo=0;

        theNucleus->StartLoop();
        G4Nucleon * pNucleon =0;

        while( (pNucleon = theNucleus->GetNextNucleon()) ) {if(NucleonNo == theCurrent) break; NucleonNo++;} 

        if ( pNucleon ) {

	  G4QGSMSplitableHadron* aTarget = new G4QGSMSplitableHadron(*pNucleon);
          pNucleon->Hit(aTarget);

	  if ( (0.06 > G4UniformRand() &&(ModelMode==SOFT)) || (ModelMode==DIFFRACTIVE ) )
	  {     // Diffractive interaction
      		G4InteractionContent * aInteraction = new G4InteractionContent(theProjectileSplitable);
      		theProjectileSplitable->SetStatus(1*theProjectileSplitable->GetStatus());

      		aInteraction->SetTarget(aTarget);
      		aInteraction->SetTargetNucleon(pNucleon);
      		aTarget->SetCollisionCount(0);
      		aTarget->SetStatus(1);                             // Mark that is Diffr. interaction

      		aInteraction->SetNumberOfDiffractiveCollisions(1);
      		aInteraction->SetNumberOfSoftCollisions(0);
      		aInteraction->SetStatus(1);

      		theInteractions.push_back(aInteraction);
		totalCuts += 1;
	  }
	  else
	  {
		// nondiffractive soft interaction occurs
		aTarget->IncrementCollisionCount(1);
	        aTarget->SetStatus(0);
        	theTargets.push_back(aTarget);

		theProjectileSplitable->IncrementCollisionCount(1);
        	theProjectileSplitable->SetStatus(0*theProjectileSplitable->GetStatus());

		G4InteractionContent * aInteraction = 
                                           new G4InteractionContent(theProjectileSplitable);
		aInteraction->SetTarget(aTarget);
        	aInteraction->SetTargetNucleon(pNucleon);
		aInteraction->SetNumberOfSoftCollisions(1);
        	aInteraction->SetStatus(0);                        // Mark that is non-Diffr. interaction
		theInteractions.push_back(aInteraction);
		totalCuts += 1;
	  }
        }
	return theProjectileSplitable;
}
