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
// $Id: G4OpAbsorption.cc 106117 2017-09-13 10:23:20Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Absorption Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        G4OpAbsorption.cc
// Description: Discrete Process -- Absorption of Optical Photons  
// Version:     1.0
// Created:     1996-05-21
// Author:      Juliet Armstrong
// Updated:     2005-07-28 - add G4ProcessType to constructor
//              2000-09-18 by Peter Gumplinger
//              > comment out warning - "No Absorption length specified" 
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
//              1998-08-25 by Stefano Magni
//              > Change process to use G4MaterialPropertiesTables
//              1998-09-03 by Peter Gumplinger
//              > Protect G4MaterialPropertyVector* AttenuationLengthVector
// mail:        gum@triumf.ca
//              magni@mi.infn.it
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpProcessSubType.hh"

#include "G4OpAbsorption.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4OpAbsorption::operator=(const G4OpAbsorption &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4OpAbsorption::G4OpAbsorption(const G4String& processName, G4ProcessType type)
              : G4VDiscreteProcess(processName, type)
{
        if (verboseLevel>0) {
           G4cout << GetProcessName() << " is created " << G4endl;
        }

        SetProcessSubType(fOpAbsorption);
}

// G4OpAbsorption::G4OpAbsorption(const G4OpAbsorpton &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4OpAbsorption::~G4OpAbsorption(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// -------------
//
G4VParticleChange*
G4OpAbsorption::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
        aParticleChange.Initialize(aTrack);

        const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        G4double thePhotonMomentum = aParticle->GetTotalMomentum();

        aParticleChange.ProposeLocalEnergyDeposit(thePhotonMomentum);

        aParticleChange.ProposeTrackStatus(fStopAndKill);

        if (verboseLevel>0) {
	   G4cout << "\n** Photon absorbed! **" << G4endl;
        }
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}


// GetMeanFreePath
// ---------------
//
G4double G4OpAbsorption::GetMeanFreePath(const G4Track& aTrack,
 				         G4double ,
				         G4ForceCondition* )
{
	const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
        const G4Material* aMaterial = aTrack.GetMaterial();

	G4double thePhotonMomentum = aParticle->GetTotalMomentum();

	G4MaterialPropertiesTable* aMaterialPropertyTable;
	G4MaterialPropertyVector* AttenuationLengthVector;
	
        G4double AttenuationLength = DBL_MAX;

	aMaterialPropertyTable = aMaterial->GetMaterialPropertiesTable();

	if ( aMaterialPropertyTable ) {
	   AttenuationLengthVector = aMaterialPropertyTable->
                                                GetProperty(kABSLENGTH);
           if ( AttenuationLengthVector ){
             AttenuationLength = AttenuationLengthVector->
                                         Value(thePhotonMomentum);
           }
           else {
//             G4cout << "No Absorption length specified" << G4endl;
           }
        } 
        else {
//           G4cout << "No Absorption length specified" << G4endl;
        }

        return AttenuationLength;
}
