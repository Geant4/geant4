// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05SpecialCuts.cc,v 1.3 1999-12-15 14:49:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//	For information related to this code contact:
//	CERN, CN Division, ASD group
// ------------------------------------------------------------
//                  15 April 1998  M.Maire
// ------------------------------------------------------------

#include "ExN05SpecialCuts.hh"
#include "G4VParticleChange.hh"
#include "G4Track.hh"
#include "G4Step.hh"

ExN05SpecialCuts::ExN05SpecialCuts(const G4String& aName)
  : G4VProcess(aName)
{
   if (verboseLevel>1) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

ExN05SpecialCuts::~ExN05SpecialCuts() 
{                                     
}                                     

G4VParticleChange* ExN05SpecialCuts::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step& 
			    )
//
// Stop the current particle, if requested by G4UserLimits 
// 			    			    			    
{
   aParticleChange.Initialize(aTrack);
   aParticleChange.SetEnergyChange(0.) ;
   aParticleChange.SetLocalEnergyDeposit (aTrack.GetKineticEnergy()) ;
   aParticleChange.SetStatusChange(fStopButAlive);
   return &aParticleChange;
}

G4double ExN05SpecialCuts::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition
                            )
{
  return DBL_MAX;
}

