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
////////////////////////////////////////////////////////////////////////
//
// File G4OpMieHG.hh
// Description: Discrete Process -- Mie Scattering of Optical Photons
// Created: 2010-07-03
// Author: Xin Qian
// Based on work from Vlasios Vasileiou
//
// This subroutine will mimic the Mie scattering based on 
// Henyey-Greenstein phase function
// Forward and backward angles are treated separately.
//
// mail: gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4OpMieHG_h
#define G4OpMieHG_h 1

#include "G4VDiscreteProcess.hh"
#include "G4OpticalPhoton.hh"

class G4OpMieHG : public G4VDiscreteProcess
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4OpMieHG(const G4String& processName = "OpMieHG",
                           G4ProcessType type = fOptical);
        ~G4OpMieHG();

private:

        G4OpMieHG(const G4OpMieHG &right);

        //////////////
        // Operators
        //////////////

        G4OpMieHG& operator=(const G4OpMieHG &right);

public:

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);  
        // Returns true -> 'is applicable' only for an optical photon.

        G4double GetMeanFreePath(const G4Track& aTrack,
                                 G4double,
                                 G4ForceCondition* );
        // Return the mean free path of Mie scattering

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                        const G4Step&  aStep);
        // This is the method implementing Mie scattering.
};
  

inline
G4bool G4OpMieHG::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
}

#endif /* G4OpMieHG_h */
