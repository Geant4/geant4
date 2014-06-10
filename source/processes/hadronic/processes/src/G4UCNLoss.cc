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
// $Id: G4UCNLoss.cc 69576 2013-05-08 13:48:13Z gcosmo $
//
////////////////////////////////////////////////////////////////////////
// Ultra Cold Neutron Loss Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:         G4UCNLoss.cc
// Description:  Discrete Process -- Loss of UCN
//               energy-independent loss cross section inside a material
// Version:      1.0
// Created:      2014-05-05
// Author:       Peter Gumplinger
// Adopted from: UCN simple absorption, Peter Fierlinger 7.9.04
//                                      Marcin Kuzniaka 21.4.06
// Updated:
// mail:         gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4UCNProcessSubType.hh"

#include "G4UCNLoss.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////

// G4UCNLoss::operator=(const G4UCNLoss &right)
// {
// }

        /////////////////
        // Constructors
        /////////////////

G4UCNLoss::G4UCNLoss(const G4String& processName, G4ProcessType type)
         : G4VDiscreteProcess(processName, type)
{
  if (verboseLevel>0) G4cout << GetProcessName() << " is created " << G4endl;

  SetProcessSubType(fUCNLoss);
}

// G4UCNLoss::G4UCNLoss(const G4UCNLoss &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4UCNLoss::~G4UCNLoss(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// -------------

G4VParticleChange*
G4UCNLoss::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  aParticleChange.ProposeTrackStatus(fStopAndKill);

  if (verboseLevel>0) G4cout << "\n** UCN lost! **" << G4endl;

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------

G4double G4UCNLoss::GetMeanFreePath(const G4Track& aTrack,
 				    G4double ,
				    G4ForceCondition* )
{
  G4double AttenuationLength = DBL_MAX;

  const G4Material* aMaterial = aTrack.GetMaterial();
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                     aMaterial->GetMaterialPropertiesTable();

  G4double crossect = 0.0;
  if (aMaterialPropertiesTable) {
     crossect = aMaterialPropertiesTable->GetConstProperty("LOSSCS");
//     if (crossect == 0.0)
//       G4cout << "No Loss Cross Section specified" << G4endl;
  }
//  else G4cout << "No Loss Cross Section specified" << G4endl;

  if (crossect) {

     // Calculate a UCN absorption length for this cross section

     G4double density = aMaterial->GetTotNbOfAtomsPerVolume();

     // Calculate cross section for a constant loss 

     crossect *= barn;
  
     // Attenuation length

     AttenuationLength = 1./density/crossect;
  }

  return AttenuationLength;
}
