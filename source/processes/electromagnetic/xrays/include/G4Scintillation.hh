// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Scintillation.hh,v 1.2 1999-10-30 00:40:01 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Scintillation Light Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:        G4Scintillation.hh  
// Description:	Discrete Process - Generation of Scintillation Photons
// Version:     1.0
// Created:     1998-11-07
// Author:      Peter Gumplinger
// Updated:     1999-10-29 add method and class decriptors
//
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4Scintillation_h
#define G4Scintillation_h 1

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
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

// Class Description:
// Discrete Process - Generation of Scintillation Photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4Scintillation : public G4VDiscreteProcess
{

private:

        //////////////
        // Operators
        //////////////

	// G4Scintillation& operator=(const G4Scintillation &right);

public: // Without description

	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

	G4Scintillation(const G4String& processName = "Scintillation");

	// G4Scintillation(const G4Scintillation &right);

	~G4Scintillation();	

        ////////////
        // Methods
        ////////////

public: // With description

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable', for any particle type.

	G4double GetMeanFreePath(const G4Track& aTrack,
				       G4double ,
                                       G4ForceCondition* );
        // Returns infinity; i. e. the process does not limit the step,
        // but sets the 'Forced' condition for the DoIt to be invoked at 
        // every step.

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
			                const G4Step&  aStep);
        // This is the method implementing the scintillation process.

	void SetTrackSecondariesFirst(const G4bool state);
        // If set, the primary particle tracking is interrupted and any
        // produced scintillation photons are tracked next. When all 
        // have been tracked, the tracking of the primary resumes.

        G4bool GetTrackSecondariesFirst() const;
        // Returns the boolean flag for tracking secondaries first.
	
        void SetScintillationYield(const G4double yield);
        // Called to set the scintillation photon yield per unit of
        // energy deposition. IMPORTANT: 'yield' has units of 1/energy
        // and is not a simple number!
        G4double GetScintillationYield() const;
        // Returns the photon yield per energy unit.

        void SetResolutionScale(const G4double scale);
        // Called to set sigma of the Gaussian photon yield distribution.
        G4double GetResolutionScale() const;
        // Returns the sigma of the Gaussian photon yield distribution.

        void SetScintillationTime(const G4double time);
        // Called to set the exponential time constant for the delayed
        // scintillation photon emission.
        G4double GetScintillationTime() const;
        // Returns the exponential time constant for the delayed 
        // scintillation photon emission.

        G4PhysicsTable* GetPhysicsTable() const;
        // Returns the address of the physics table.

        void DumpPhysicsTable() const;
        // Prints the physics table.

private:

        void BuildThePhysicsTable();

        ///////////////////////
        // Class Data Members
        ///////////////////////

protected:

        G4PhysicsTable* thePhysicsTable;
        //  A Physics Table can be either a cross-sections table or
        //  an energy table (or can be used for other specific
        //  purposes).

private:

	G4bool fTrackSecondariesFirst;

        G4double ScintillationYield;
        G4double ScintillationTime;
        G4double ResolutionScale;

};

////////////////////
// Inline methods
////////////////////

inline 
G4bool G4Scintillation::IsApplicable(const G4ParticleDefinition& aParticleType)
{
        return true;
}

inline 
void G4Scintillation::SetTrackSecondariesFirst(const G4bool state) 
{ 
	fTrackSecondariesFirst = state;
}

inline
G4bool G4Scintillation::GetTrackSecondariesFirst() const
{
        return fTrackSecondariesFirst;
}

inline
void G4Scintillation::SetScintillationYield(const G4double yield)
{
        ScintillationYield = yield;
}

inline
G4double G4Scintillation::GetScintillationYield() const
{
        return ScintillationYield;
}

inline
void G4Scintillation::SetResolutionScale(const G4double scale)
{
        ResolutionScale = scale;
}

inline
G4double G4Scintillation::GetResolutionScale() const
{
        return ResolutionScale;
}

inline
void G4Scintillation::SetScintillationTime(const G4double time)
{
        ScintillationTime = time;
}

inline
G4double G4Scintillation::GetScintillationTime() const
{
        return ScintillationTime;
}

inline
G4PhysicsTable* G4Scintillation::GetPhysicsTable() const
{
        return thePhysicsTable;
}

inline
void G4Scintillation::DumpPhysicsTable() const
{
        G4int PhysicsTableSize = thePhysicsTable->entries();
        G4PhysicsOrderedFreeVector *v;

        for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
        {
        	v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
        	v->DumpValues();
        }
}

#endif /* G4Scintillation_h */
