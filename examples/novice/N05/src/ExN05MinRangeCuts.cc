// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05MinRangeCuts.cc,v 1.1 1999-01-07 16:06:17 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD Group
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
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
     G4cout << GetProcessName() << " is created "<< endl;
   }
   SetProcessType(fUserDefined);
}

ExN05MinRangeCuts::~ExN05MinRangeCuts()
{}

ExN05MinRangeCuts::ExN05MinRangeCuts(ExN05MinRangeCuts& right)
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
     G4Material* aMaterial = aTrack.GetMaterial();
     G4double    rangeNow = DBL_MAX;
     if (aParticleDef->GetPDGCharge() != 0.0) {
       rangeNow = G4EnergyLossTables::GetRange(aParticleDef,eKine,aMaterial);
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
