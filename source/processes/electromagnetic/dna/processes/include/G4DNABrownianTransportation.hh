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
// Author: Mathieu Karamitros, kara@cenbg.in2p3.fr

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#ifndef G4ITBROWNIANTRANSPORTATION_H
#define G4ITBROWNIANTRANSPORTATION_H

#include "G4ITTransportation.hh"

class G4SafetyHelper;

/* \brief { The transportation method implemented is the one from
 *         Ermak-McCammon : J. Chem. Phys. 69, 1352 (1978)}
 */

class G4DNABrownianTransportation : public G4ITTransportation
{
public:
	G4DNABrownianTransportation(const G4String& aName =  "DNABrownianTransportation",
	                            G4int verbosityLevel= 0);
	G4IT_ADD_CLONE(G4VITProcess,G4DNABrownianTransportation)
	virtual ~G4DNABrownianTransportation();
	G4DNABrownianTransportation(const G4DNABrownianTransportation& other);
	G4DNABrownianTransportation& operator=(const G4DNABrownianTransportation& other);

	virtual void BuildPhysicsTable(const G4ParticleDefinition&);

	virtual void StartTracking(G4Track* aTrack);

	virtual void ComputeStep(const G4Track&,
			const G4Step&,
			const double,
			double&) ;

	virtual G4double AlongStepGetPhysicalInteractionLength( const G4Track& /*track*/,
			G4double /*previousStepSize*/,
			G4double /*currentMinimumStep*/,
			G4double& /*currentSafety*/,
			G4GPILSelection* /*selection*/);
	virtual G4VParticleChange* PostStepDoIt( const G4Track& track, const G4Step& ) ;

	virtual G4VParticleChange* AlongStepDoIt(const G4Track& track, const G4Step&);

protected:
	void Diffusion(const G4Track& track);

	//________________________________________________________________
	// Process information
	struct G4ITBrownianState : public G4ITTransportationState
	{
	public :
		G4ITBrownianState();
		virtual ~G4ITBrownianState(){;}
		virtual G4String GetType()
		{
			return "G4ITBrownianState";
		}

		G4bool  fPathLengthWasCorrected;
		G4bool fTimeStepReachedLimit;
	};


	G4bool fUseMaximumTimeBeforeReachingBoundary;
	G4Material* fNistWater ;

	// Water density table
	const std::vector<G4double>* fpWaterDensity;
};

#endif // G4ITBROWNIANTRANSPORTATION_H
