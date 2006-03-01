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
// $Id: G4StepLimiterPerRegion.cc,v 1.2 2006-03-01 13:52:01 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4StepLimiterPerRegion.hh"
// #include "G4StepLimiterMessenger.hh"
#include "G4VPhysicalVolume.hh"
// #include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StepLimiterPerRegion::G4StepLimiterPerRegion(const G4String& processName)
 : G4VDiscreteProcess(processName),
   MaxChargedStep(DBL_MAX)
{
  //   pMess = new G4StepLimiterMessenger(this);
  //   gasVolume = Histo::GetPointer()->GasVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4StepLimiterPerRegion::~G4StepLimiterPerRegion() 
{ 
  // delete pMess; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4StepLimiterPerRegion::IsApplicable(const G4ParticleDefinition& particle)
{
  return (particle.GetPDGCharge() != 0. && !(particle.IsShortLived()));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4StepLimiterPerRegion::SetMaxStep(G4double step) 
{
  MaxChargedStep = step;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4StepLimiterPerRegion::PostStepGetPhysicalInteractionLength(
                                              const G4Track& aTrack,
                                                    G4double,
                                                    G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  if(aTrack.GetVolume()->GetName() == "Mylar" ) ProposedStep = 1.e-7*mm;
  else                                          ProposedStep = MaxChargedStep;

  return ProposedStep;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VParticleChange* G4StepLimiterPerRegion::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
  aParticleChange.Initialize(aTrack);

  if(aTrack.GetVolume()->GetName() == "Mylar" ) // "Pipe1" ) 
  {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
  }
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
