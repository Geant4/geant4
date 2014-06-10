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
// $Id: G4UCNAbsorption.hh 71487 2013-06-17 08:19:40Z gcosmo $
//
//
////////////////////////////////////////////////////////////////////////
// Ultra Cold Neutron (UCN) Absortion Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:         G4UCNAbsorption.hh
// Description:  Discrete Process -- Bulk absorption of UCN
// Version:      1.0
// Created:      2014-05-05
// Author:       Peter Gumplinger
// Adopted from: UCNSimpleAbsorption by Peter Fierlinger 7.9.2004
// Updated:
// mail:         gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef G4UCNABSORPTION_HH
#define G4UCNABSORPTION_HH 1

/////////////
// Includes
/////////////

#include "G4VDiscreteProcess.hh"

#include "G4Neutron.hh"

// Class Description:
// Discrete Process -- Bulk loss of Ultra Cold Neutrons.
// Absorbs UCN due to absorption cross section of the material.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////

class G4UCNAbsorption : public G4VDiscreteProcess
{

public:

        ////////////////////////////////
        // Constructors and Destructor
        ////////////////////////////////

        G4UCNAbsorption(const G4String& processName = "UCNAbsorption",
                                 G4ProcessType type = fUCN);
	virtual ~G4UCNAbsorption();

private:

        G4UCNAbsorption(const G4UCNAbsorption &right);

        //////////////
        // Operators
        //////////////

        G4UCNAbsorption& operator=(const G4UCNAbsorption &right);

public:

        ////////////
        // Methods
        ////////////

        G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
        // Returns true -> 'is applicable' only for an UCN.

	G4double GetMeanFreePath(const G4Track& aTrack,
				 G4double ,
				 G4ForceCondition* condition);

	G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
 				        const G4Step&  aStep);

private:

};

////////////////////
// Inline methods
////////////////////

inline
G4bool G4UCNAbsorption::IsApplicable(const G4ParticleDefinition& aParticleType)
{
   return ( &aParticleType == G4Neutron::NeutronDefinition() );
}

#endif /* G4UCNABSORPTION_HH */
