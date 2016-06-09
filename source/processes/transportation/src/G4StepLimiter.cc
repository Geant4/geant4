//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StepLimiter.cc,v 1.1 2004/07/26 00:42:59 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// --------------------------------------------------------------
// History
//
// 23-01-04 first implementation  (H.Kurashige)
// --------------------------------------------------------------

#include "G4StepLimiter.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"


////////////////////////////////////
G4StepLimiter::G4StepLimiter(const G4String& aName)
  : G4VProcess(aName)
{
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













