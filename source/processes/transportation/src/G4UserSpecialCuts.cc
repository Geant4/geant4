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
//
// --------------------------------------------------------------
// History
//
// 15-04-98 first implementation, mma
// 07-04-03 migrade to cut per region (V.Ivanchenko)
// 18-09-03 substitute manager for the loss tables (V.Ivanchenko)
// 23-01-04 add protection for charged geantino in range cut (H.Kurashige)
// 09-09-04 tracking cut applied only if Rmin or Emin > DBL_MIN
// 21-01-11 changed order of checks: 1st energy tracking cut, 2nd track
//          length, 3d time cut, 4th range for charged particles 
//          with non-zero mass; removed string comparisons (V.Ivanchenko)
// --------------------------------------------------------------

#include "G4UserSpecialCuts.hh"
#include "G4TransportationProcessType.hh"

#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4UserLimits.hh"
#include "G4VParticleChange.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UserSpecialCuts::G4UserSpecialCuts(const G4String& aName)
  : G4VProcess(aName, fGeneral  )
{
  // set Process Sub Type
  SetProcessSubType(static_cast<int>(USER_SPECIAL_CUTS));

  if (verboseLevel>0)
    {
      G4cout << GetProcessName() << " is created "<< G4endl;
    }
  theLossTableManager = G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UserSpecialCuts::~G4UserSpecialCuts()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4UserSpecialCuts::G4UserSpecialCuts(G4UserSpecialCuts& right)
  : G4VProcess(right)
{
  theLossTableManager = G4LossTableManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4double G4UserSpecialCuts::
PostStepGetPhysicalInteractionLength( const G4Track& aTrack,
                                            G4double, // previousStepSize
                                            G4ForceCondition* condition  )
{
  // condition is set to "Not Forced"
  *condition = NotForced;

   G4double ProposedStep = DBL_MAX;
   G4UserLimits* pUserLimits =
                 aTrack.GetVolume()->GetLogicalVolume()->GetUserLimits();
   if (pUserLimits)
   {
     // check max kinetic energy first
     //
     G4double Ekine = aTrack.GetKineticEnergy();
     if(Ekine <= pUserLimits->GetUserMinEkine(aTrack)) { return 0.; }

     // max track length
     //
     ProposedStep = (pUserLimits->GetUserMaxTrackLength(aTrack)
                   - aTrack.GetTrackLength());
     if (ProposedStep < 0.) { return 0.; }

     // max time limit
     //
     G4double tlimit = pUserLimits->GetUserMaxTime(aTrack);
     if(tlimit < DBL_MAX) {
       G4double beta  = (aTrack.GetDynamicParticle()->GetTotalMomentum())
	 /(aTrack.GetTotalEnergy());
       G4double dTime = (tlimit - aTrack.GetGlobalTime());
       G4double temp  = beta*c_light*dTime;
       if (temp < 0.) { return 0.; }
       if (ProposedStep > temp) { ProposedStep = temp; }
     }
                 
     // min remaining range 
     // (only for charged particle except for chargedGeantino)
     //
     G4double Rmin = pUserLimits->GetUserMinRange(aTrack);
     if (Rmin > DBL_MIN) {
       G4ParticleDefinition* Particle = aTrack.GetDefinition();
       if ( (Particle->GetPDGCharge() != 0.) && (Particle->GetPDGMass() > 0.0))
       {
	 const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
	 G4double RangeNow = theLossTableManager->GetRange(Particle,Ekine,couple);
         G4double temp = RangeNow - Rmin;
         if (temp < 0.) { return 0.; }
         if (ProposedStep > temp) { ProposedStep = temp; }
       }	 
     }
   }
   return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange*
G4UserSpecialCuts::PostStepDoIt( const G4Track& aTrack,
                                 const G4Step&  )
//
// Kill the current particle, if requested by G4UserLimits
//
{
   aParticleChange.Initialize(aTrack);
   aParticleChange.ProposeEnergy(0.) ;
   aParticleChange.ProposeLocalEnergyDeposit(aTrack.GetKineticEnergy()) ;
   aParticleChange.ProposeTrackStatus(fStopAndKill);
   return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......












