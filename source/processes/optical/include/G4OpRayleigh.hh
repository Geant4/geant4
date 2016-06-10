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
// $Id: G4OpRayleigh.hh 84596 2014-10-17 07:34:41Z gcosmo $
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
// Updated:     2014-08-20 allow for more material types
//              2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1997-04-09 by Peter Gumplinger
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

// Class Description:
// Discrete Process -- Rayleigh scattering of optical photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4OpRayleigh : public G4VDiscreteProcess 
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////
 
        G4OpRayleigh(const G4String& processName = "OpRayleigh",
                              G4ProcessType type = fOptical);
	~G4OpRayleigh();

private:

        G4OpRayleigh(const G4OpRayleigh &right);

        //////////////
        // Operators
        //////////////

        G4OpRayleigh& operator=(const G4OpRayleigh &right);

public:

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

        void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
        // Build thePhysicsTable at a right time

        G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double ,
                                 G4ForceCondition* );
        // Returns the mean free path for Rayleigh scattering in water.
        // --- Not yet implemented for other materials! ---

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                       const G4Step&  aStep);
        // This is the method implementing Rayleigh scattering.

        G4PhysicsTable* GetPhysicsTable() const;
        // Returns the address of the physics table.

        void DumpPhysicsTable() const;
        // Prints the physics table.

private:

        /////////////////////
        // Helper Functions
        /////////////////////

        /// Calculates the mean free paths for a material as a function of 
        /// photon energy
        ///
        /// @param[in] material information
        /// @return the mean free path vector
        G4PhysicsOrderedFreeVector* 
        CalculateRayleighMeanFreePaths( const G4Material* material ) const;

        ///////////////////////
        // Class Data Members
        ///////////////////////

protected:

        G4PhysicsTable* thePhysicsTable;
        //  A Physics Table can be either a cross-sections table or
        //  an energy table (or can be used for other specific
        //  purposes).

private:
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
