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
// $Id: ExN05MinRangeCuts.cc,v 1.7 2006/06/29 17:53:28 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file
// --------------------------------------------------------------
//                   15 April 1998 M.Maire
// --------------------------------------------------------------

#include "ExN05MinRangeCuts.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

ExN05MinRangeCuts::ExN05MinRangeCuts(const G4String& aName)
  : ExN05SpecialCuts(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
   SetProcessType(fUserDefined);
}

ExN05MinRangeCuts::~ExN05MinRangeCuts()
{}

ExN05MinRangeCuts::ExN05MinRangeCuts(ExN05MinRangeCuts&)
  : ExN05SpecialCuts()
{}


G4double ExN05MinRangeCuts::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
			     G4double ,
			     G4ForceCondition* condition
			    )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double     proposedStep = DBL_MAX;
   // get the pointer to UserLimits
   G4UserLimits* pUserLimits = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
   const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
   G4ParticleDefinition* aParticleDef = aTrack.GetDefinition();

   if (pUserLimits) {
     G4double temp = DBL_MAX;

     //min remaining range
     G4double    eKine     = aParticle->GetKineticEnergy();
     const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
     G4double    rangeNow = DBL_MAX;
     if (aParticleDef->GetPDGCharge() != 0.0) {
       rangeNow = G4EnergyLossTables::GetRange(aParticleDef,eKine,couple);
     }
     temp = (rangeNow - pUserLimits->GetUserMinRange(aTrack));
     if (temp < 0.) {
       proposedStep = 0.;
     } else {
       if (proposedStep > temp) proposedStep = temp;
     }

   }
   return proposedStep;
}
