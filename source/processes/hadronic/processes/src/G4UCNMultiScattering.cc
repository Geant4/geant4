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
// $Id: G4UCNMultiScattering.cc 69576 2013-05-08 13:48:13Z gcosmo $
//
///////////////////////////////////////////////////////////////////////
// UCN Multiple Scattering Class Implementation
///////////////////////////////////////////////////////////////////////
//
// File:        G4UCNMultiScattering.cc
// Description: G4VDiscreteProcess -- MultiScattering of UCNs
// Version:     1.0
// Created:     2014-05-12
// Author:      Peter Gumplinger
//              adopted from Geant4UCN by Peter Fierlinger (7.9.04) and
//              Marcin Kuzniak (21.4.06)
//              Calculate "elastic scattering" inside materials
// Updated:
//
// mail:        gum@triumf.ca
//
///////////////////////////////////////////////////////////////////////

#include "G4UCNProcessSubType.hh"

#include "G4UCNMultiScattering.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

/////////////////////////
// Class Implementation
/////////////////////////

        //////////////
        // Operators
        //////////////
 
        // G4UCNMultiScattering::operator=(const G4UCNMultiScattering &right)
        // {
        // }

        /////////////////
        // Constructors
        /////////////////

G4UCNMultiScattering::G4UCNMultiScattering(const G4String& processName,
                                           G4ProcessType type)
                    : G4VDiscreteProcess(processName, type)
{
  if (verboseLevel>0) G4cout << GetProcessName() << " is created " << G4endl;

  SetProcessSubType(fUCNMultiScattering);
}

// G4UCNMultiScattering::G4UCNMultiScattering(const G4UCNMultiScattering &right)
// {
// }

        ////////////////
        // Destructors
        ////////////////

G4UCNMultiScattering::~G4UCNMultiScattering(){}

        ////////////
        // Methods
        ////////////

// PostStepDoIt
// ------------

G4VParticleChange*
G4UCNMultiScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  if ( verboseLevel > 0 ) G4cout << "UCNMULTISCATTER at: "
     << aTrack.GetProperTime()/s << "s, "
     << aTrack.GetGlobalTime()/s << "s. "
     << ", after track length " << aTrack.GetTrackLength()/cm << "cm, "
           << "in volume "
           << aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()
           << G4endl;

  G4ThreeVector scattered = Scatter();
      
  aParticleChange.ProposeMomentumDirection(-scattered);

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

// GetMeanFreePath
// ---------------

G4double G4UCNMultiScattering::GetMeanFreePath(const G4Track& aTrack,
                                               G4double ,
                                               G4ForceCondition*)
{
  G4double AttenuationLength = DBL_MAX;

  const G4Material* aMaterial = aTrack.GetMaterial();
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                                     aMaterial->GetMaterialPropertiesTable();

  G4double crossect = 0.0;
  if (aMaterialPropertiesTable) {
     crossect = aMaterialPropertiesTable->GetConstProperty("SCATCS");
//     if(crossect == 0.0)
//       G4cout << "No UCN MultiScattering length specified" << G4endl;
  }
//  else G4cout << "No UCN MultiScattering length specified" << G4endl;

  if (crossect) {

    // Calculate a UCN MultiScattering length for this cross section

    G4double density = aMaterial->GetTotNbOfAtomsPerVolume();

    crossect *= barn;
  
    // attenuation length in mm
    AttenuationLength = 1./density/crossect; 
  }

  return AttenuationLength;
}

G4ThreeVector G4UCNMultiScattering::Scatter()
{

  G4ThreeVector final(0.,0.,1.);

  // Make a simple uniform distribution in 4 pi
  // apply scattering, calculate angle phi, theta

  G4double theta = std::acos(2*G4UniformRand()-1);
  G4double phi = G4UniformRand() * 2 * pi;

  final.rotateY(theta);
  final.rotateZ(phi);
  final = final.unit();

  return final;
}
