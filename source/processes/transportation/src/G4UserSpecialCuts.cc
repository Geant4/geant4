// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UserSpecialCuts.cc,v 1.3 1999-12-15 14:53:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// --------------------------------------------------------------
// History
//
// 15-04-98 first implementation, mma                   
// --------------------------------------------------------------

#include "G4UserSpecialCuts.hh"

#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UserSpecialCuts::G4UserSpecialCuts(const G4String& aName)
  : G4VProcess(aName)
{
   if (verboseLevel>0) {
     G4cout << GetProcessName() << " is created "<< G4endl;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UserSpecialCuts::~G4UserSpecialCuts()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UserSpecialCuts::G4UserSpecialCuts(G4UserSpecialCuts& right)
  : G4VProcess(right)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
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
       //min remaining range (only for charged particle)
       G4ParticleDefinition* Particle = aTrack.GetDefinition();
       if (Particle->GetPDGCharge() != 0.)
         {G4double              Ekine    = aTrack.GetKineticEnergy();
          G4Material*           Material = aTrack.GetMaterial();
          G4double RangeNow = G4EnergyLossTables::GetRange(Particle,Ekine,Material);
          temp = (RangeNow - pUserLimits->GetUserMinRange(aTrack));
          if (temp < 0.) return 0.;
          if (ProposedStep > temp) ProposedStep = temp;
          //min kinetic energy (only for charged particle)
          G4double Emin = pUserLimits->GetUserMinEkine(aTrack);
          G4double Rmin = G4EnergyLossTables::GetRange(Particle,Emin,Material);
          temp = RangeNow - Rmin;
          if (temp < 0.) return 0.;
          if (ProposedStep > temp) ProposedStep = temp;
	 }         
     }   
   return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
