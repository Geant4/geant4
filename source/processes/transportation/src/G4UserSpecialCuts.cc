// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserSpecialCuts.cc,v 1.1 1999-01-07 16:14:11 gunter Exp $
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

#include "G4UserSpecialCuts.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

G4UserSpecialCuts::G4UserSpecialCuts(const G4String& aName)
  : G4VProcess(aName)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< endl;
   }
}

G4UserSpecialCuts::~G4UserSpecialCuts()
{}

G4UserSpecialCuts::G4UserSpecialCuts(G4UserSpecialCuts& right)
  : G4VProcess(right)
{}

 
G4double G4UserSpecialCuts::PostStepGetPhysicalInteractionLength(
                             const G4Track& aTrack,
			     G4double   previousStepSize,
			     G4ForceCondition* condition
			    )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double ProposedStep = DBL_MAX;
   G4UserLimits* pUserLimits = aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
   if (pUserLimits)
     { //max track length
       ProposedStep = (pUserLimits->GetUserMaxTrackLength(aTrack) - aTrack.GetTrackLength());
       if (ProposedStep < 0.) return 0.;
       //max time limit
       G4double beta = (aTrack.GetDynamicParticle()->GetTotalMomentum())/(aTrack.GetTotalEnergy());
       G4double dTime= (pUserLimits->GetUserMaxTime(aTrack) - aTrack.GetGlobalTime());
       G4double temp = beta*c_light*dTime;
       if (temp < 0.) return 0.;
       if (ProposedStep > temp) ProposedStep = temp;                  
       //min remaining range
       G4ParticleDefinition* Particle = aTrack.GetDefinition();
       G4double              Ekine    = aTrack.GetKineticEnergy();
       G4Material*           Material = aTrack.GetMaterial();
       G4double RangeNow = G4EnergyLossTables::GetRange(Particle,Ekine,Material);
       temp = (RangeNow - pUserLimits->GetUserMinRange(aTrack));
       if (temp < 0.) return 0.;
       if (ProposedStep > temp) ProposedStep = temp;
       //min kinetic energy
       G4double Emin = pUserLimits->GetUserMinEkine(aTrack);
       G4double Rmin = G4EnergyLossTables::GetRange(Particle,Emin,Material);
       temp = RangeNow - Rmin;
       if (temp < 0.) return 0.;
       if (ProposedStep > temp) ProposedStep = temp;        
     }   
   return ProposedStep;
}

G4VParticleChange* G4UserSpecialCuts::PostStepDoIt(
			     const G4Track& aTrack,
			     const G4Step& 
			    )
//
// Kill the current particle, if requested by G4UserLimits 
// 			    			    			    
{
   aParticleChange.Initialize(aTrack);
   aParticleChange.SetEnergyChange(0.) ;
   aParticleChange.SetLocalEnergyDeposit (aTrack.GetKineticEnergy()) ;
   aParticleChange.SetStatusChange(fStopAndKill);
   return &aParticleChange;
}
