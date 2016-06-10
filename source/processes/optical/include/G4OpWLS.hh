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
// $Id: G4OpWLS.hh 85354 2014-10-28 09:58:12Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon WaveLength Shifting (WLS) Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpWLS.hh
// Description: Discrete Process -- Wavelength Shifting of Optical Photons 
// Version:     1.0
// Created:     2003-05-13
// Author:      John Paul Archambault
//              (Adaptation of G4Scintillation and G4OpAbsorption)
// Updated:     2005-07-28 add G4ProcessType to constructor
//              2006-05-07 - add G4VWLSTimeGeneratorProfile
// mail:        gum@triumf.ca
//              jparcham@phys.ualberta.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpWLS_h
#define G4OpWLS_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include "G4VWLSTimeGeneratorProfile.hh"

// Class Description:
// Discrete Process -- Bulk absorption of Optical Photons.
// Class inherits publicly from G4VDiscreteProcess
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4VWLSTimeGeneratorProfile;

class G4OpWLS : public G4VDiscreteProcess 
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4OpWLS(const G4String& processName = "OpWLS",
                         G4ProcessType type = fOptical);
        ~G4OpWLS();

private:

        G4OpWLS(const G4OpWLS &right);

        //////////////
        // Operators
        //////////////

        G4OpWLS& operator=(const G4OpWLS &right);

public:

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an optical photon.

        void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
        // Build the WLS integral table at the right time

        G4double GetMeanFreePath(const G4Track& aTrack,
                                 G4double ,
                                 G4ForceCondition* );
        // Returns the absorption length for bulk absorption of optical
        // photons in media with a specified attenuation length.

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                        const G4Step&  aStep);
        // This is the method implementing bulk absorption of optical
        // photons.

        G4PhysicsTable* GetIntegralTable() const;
        // Returns the address of the WLS integral table.

        void DumpPhysicsTable() const;
        // Prints the WLS integral table.

        void UseTimeProfile(const G4String name);
        // Selects the time profile generator

protected:

        G4VWLSTimeGeneratorProfile* WLSTimeGeneratorProfile;
        G4PhysicsTable* theIntegralTable;

};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4OpWLS::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

inline
G4PhysicsTable* G4OpWLS::GetIntegralTable() const
{
  return theIntegralTable;
}

inline
void G4OpWLS::DumpPhysicsTable() const
{
  G4int PhysicsTableSize = theIntegralTable->entries();
  G4PhysicsOrderedFreeVector *v;
 
  for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
    {
      v = (G4PhysicsOrderedFreeVector*)(*theIntegralTable)[i];
      v->DumpValues();
    }
}

#endif /* G4OpWLS_h */
