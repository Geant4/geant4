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
// $Id: ExN05MinEkineCuts.cc,v 1.8 2006/06/29 17:53:25 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
// --------------------------------------------------------------
//                   15 April 1998 M.Maire
// --------------------------------------------------------------

#include "ExN05MinEkineCuts.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

ExN05MinEkineCuts::ExN05MinEkineCuts(const G4String& aName)
  : ExN05SpecialCuts(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
   SetProcessType(fUserDefined);
}

ExN05MinEkineCuts::~ExN05MinEkineCuts()
{}

ExN05MinEkineCuts::ExN05MinEkineCuts(ExN05MinEkineCuts&)
  : ExN05SpecialCuts()
{}

 
G4double ExN05MinEkineCuts::PostStepGetPhysicalInteractionLength(
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
     //min kinetic energy
     G4double temp = DBL_MAX;
     G4double    eKine     = aParticle->GetKineticEnergy();
     const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
     G4double eMin = pUserLimits->GetUserMinEkine(aTrack);

     G4double    rangeNow = DBL_MAX;
     if (aParticleDef->GetPDGCharge() != 0.0) {
       rangeNow = G4EnergyLossTables::GetRange(aParticleDef,eKine,couple);
     }
     if (eKine < eMin ) {
       proposedStep = 0.;
     } else if (aParticleDef->GetPDGCharge() != 0.0) {
       // charged particles only
       G4double rangeMin = G4EnergyLossTables::GetRange(aParticleDef,eMin,couple);
       temp = rangeNow - rangeMin;
       if (temp<0.) temp=0.;
       if (proposedStep > temp) proposedStep = temp;
     }
   }
   return proposedStep;
}
