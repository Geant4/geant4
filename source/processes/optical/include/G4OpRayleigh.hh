// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpRayleigh.hh,v 1.1 1999-01-07 16:14:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Optical Photon Rayleigh Scattering Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpRayleigh.hh
// Description: Discrete Process -- Rayleigh scattering of optical photons 
// Version:     1.0
// Created:     1996-05-31
// Author:      Juliet Armstrong
// Updated:     1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpRayleigh_h
#define G4OpRayleigh_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

/////////////////////
// Class Definition
/////////////////////

class G4OpRayleigh : public G4VDiscreteProcess 
{

private:
 
        //////////////
        // Operators
        //////////////

        // G4OpRayleigh& operator=(const G4OpRayleigh &right);

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////
 
        G4OpRayleigh(const G4String& processName = "Rayleigh Scattering");

        // G4OpRayleigh(const G4OpRayleigh &right);

	~G4OpRayleigh();

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

        G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double ,
                                 G4ForceCondition* );

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                       const G4Step&  aStep);

        G4PhysicsTable* GetPhysicsTable() const;
        //  Returns the address of the physics table.

        void DumpPhysicsTable() const;

private:

        void BuildThePhysicsTable();

        /////////////////////
        // Helper Functions
        /////////////////////

	G4PhysicsOrderedFreeVector* RayleighAttenuationLengthGenerator(
					G4MaterialPropertiesTable *aMPT);

        ///////////////////////
        // Class Data Members
        ///////////////////////

protected:

        G4PhysicsTable* thePhysicsTable;
        //  A Physics Table can be either a cross-sections table or
        //  an energy table (or can be used for other specific
        //  purposes).

};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4OpRayleigh::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
void G4OpRayleigh::DumpPhysicsTable() const

{
        G4int PhysicsTableSize = thePhysicsTable->entries();
        G4PhysicsOrderedFreeVector *v;

        for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
        {
                v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
                v->DumpValues();
        }
}

inline G4PhysicsTable* G4OpRayleigh::GetPhysicsTable() const
{
  return thePhysicsTable;
}


#endif /* G4OpRayleigh_h */
