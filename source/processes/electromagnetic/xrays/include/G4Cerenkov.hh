// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Cerenkov.hh,v 1.1 1999-01-07 16:11:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
////////////////////////////////////////////////////////////////////////
// Cerenkov Radiation Class Definition 
////////////////////////////////////////////////////////////////////////
//
// File:        G4Cerenkov.hh  
// Description:	Continuous Process -- Generation of Cerenkov Photons
// Version:     2.0
// Created:     1996-02-21
// Author:      Juliet Armstrong
// Updated:     1997-04-09 by Peter Gumplinger
//              > G4MaterialPropertiesTable; new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4Cerenkov_h
#define G4Cerenkov_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VContinuousProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh" 
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

/////////////////////
// Class Definition
/////////////////////

class G4Cerenkov : public G4VContinuousProcess  
{

private:

        //////////////
        // Operators
        //////////////

	// G4Cerenkov& operator=(const G4Cerenkov &right);

public:

	////////////////////////////////
	// Constructors and Destructor
	////////////////////////////////

	G4Cerenkov(const G4String& processName = "Cerenkov");

	// G4Cerenkov(const G4Cerenkov &right);

	~G4Cerenkov();	

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);

	G4double GetContinuousStepLimit(const G4Track& aTrack,
					G4double  ,
					G4double  ,
                                        G4double& );

	G4VParticleChange* AlongStepDoIt(const G4Track& aTrack, 
					const G4Step&  aStep);

	void SetTrackSecondariesFirst(const G4bool state);
	
	void SetMaxNumPhotonsPerStep(const G4int NumPhotons); 

        G4PhysicsTable* GetPhysicsTable() const;
        //  Returns the address of the physics table.

        void DumpPhysicsTable() const;

private:

        void BuildThePhysicsTable();

	/////////////////////
	// Helper Functions
	/////////////////////

	G4double GetAverageNumberOfPhotons(const G4DynamicParticle *aParticle,
		    		const G4Material *aMaterial,
				const G4MaterialPropertyVector* Rindex) const;

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
	G4int  fMaxPhotons;
};

////////////////////
// Inline methods
////////////////////

inline 
G4bool G4Cerenkov::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return (aParticleType.GetPDGCharge() != 0);
}

inline 
void G4Cerenkov::SetTrackSecondariesFirst(const G4bool state) 
{ 
	fTrackSecondariesFirst = state;
}

inline
void G4Cerenkov::SetMaxNumPhotonsPerStep(const G4int NumPhotons) 
{ 
	fMaxPhotons = NumPhotons;
}

inline
void G4Cerenkov::DumpPhysicsTable() const
{
        G4int PhysicsTableSize = thePhysicsTable->entries();
        G4PhysicsOrderedFreeVector *v;

        for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
        {
        	v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
        	v->DumpValues();
        }
}

inline G4PhysicsTable* G4Cerenkov::GetPhysicsTable() const
{
  return thePhysicsTable;
}

#endif /* G4Cerenkov_h */
