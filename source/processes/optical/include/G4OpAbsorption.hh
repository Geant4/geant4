// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpAbsorption.hh,v 1.3 1999-10-30 00:37:24 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Absorption Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpAbsorption.hh
// Description: Discrete Process -- Bulk absorption of Optical Photons
// Version:     1.0
// Created:     1996-05-21
// Author:      Juliet Armstrong
// Updated:     1999-10-29 add method and class decriptors
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
//              1998-08-25 by Stefano Magni
//              > Change process to use G4MaterialPropertiesTables
// mail:        gum@triumf.ca
//              magni@mi.infn.it
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpAbsorption_h
#define G4OpAbsorption_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"

// Class Description:
// Discrete Process -- Bulk absorption of Optical Photons.
// Class inherits publicly from G4VDiscreteProcess
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4OpAbsorption : public G4VDiscreteProcess 
{

private:

        //////////////
        // Operators
        //////////////

        // G4OpAbsorption& operator=(const G4OpAbsorption &right);

public: // Without description

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4OpAbsorption(const G4String& processName = "Absorption");

        // G4OpAbsorption(const G4OpAbsorption &right);

	~G4OpAbsorption();

	////////////
	// Methods
        ////////////

public: // With description

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

	G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double ,
				 G4ForceCondition* );
        // Returns the absorption length for bulk absorption of optical
        // photons in media with a specified attenuation length. 

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
 				        const G4Step&  aStep);
        // This is the method implementing bulk absorption of optical 
        // photons.

};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4OpAbsorption::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

#endif /* G4OpAbsorption_h */
