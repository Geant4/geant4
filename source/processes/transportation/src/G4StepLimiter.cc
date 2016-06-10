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
// $Id: G4StepLimiter.cc 68048 2013-03-13 14:34:07Z gcosmo $
//
// --------------------------------------------------------------
// History
//
// 23-01-04 first implementation  (H.Kurashige)
// --------------------------------------------------------------

#include "G4StepLimiter.hh"
#include "G4TransportationProcessType.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"


////////////////////////////////////
G4StepLimiter::G4StepLimiter(const G4String& aName)
  : G4VProcess(aName, fGeneral )
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(STEP_LIMITER));

  if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}


////////////////////////
G4StepLimiter::~G4StepLimiter()
{
}


////////////////////////
G4StepLimiter::G4StepLimiter(G4StepLimiter& right)
  : G4VProcess(right)
{
}

 
////////////////
G4double 
  G4StepLimiter::PostStepGetPhysicalInteractionLength( 
				       const G4Track& aTrack,
				       G4double, // previousStepSize
				       G4ForceCondition* condition  )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double proposedStep = DBL_MAX;
   G4UserLimits* pUserLimits =
                 aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
   if (pUserLimits) {
     // max allowed step length
     //
     proposedStep = pUserLimits->GetMaxAllowedStep(aTrack);
     if (proposedStep < 0.) proposedStep =0.; 
   }
   return proposedStep;
}

///////////////
G4VParticleChange*
  G4StepLimiter::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
// Do Nothing
//
{
   aParticleChange.Initialize(aTrack);
   return &aParticleChange;
}













