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
/// \file electromagnetic/TestEm10/src/StepMax.cc
/// \brief Implementation of the StepMax class
//
// $Id$
//
///////////////////////////////////////////////////////////////////////////

#include "StepMax.hh"
// #include "StepMaxMessenger.hh"
// #include "Histo.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

///////////////////////////////////////////////////////////////////////////

StepMax::StepMax(const G4String& processName)
 : G4VDiscreteProcess(processName),
   MaxChargedStep(DBL_MAX),
   thDensity(0.1*gram/cm3),
   first(true)
{
  //  pMess = new StepMaxMessenger(this);
  //  histo = Histo::GetPointer();
}

//////////////////////////////////////////////////////////////////////////////

StepMax::~StepMax() 
{ 
  // delete pMess; 
}

////////////////////////////////////////////////////////////////////////////

G4bool StepMax::IsApplicable(const G4ParticleDefinition& particle)
{
  return (  particle.GetPDGCharge() != 0. );
}

/////////////////////////////////////////////////////////////////////////////

void StepMax::SetMaxStep(G4double step) {MaxChargedStep = step;}

/////////////////////////////////////////////////////////////////////////////

G4double StepMax::PostStepGetPhysicalInteractionLength(
                                              const G4Track& aTrack,
                                                    G4double,
                                                    G4ForceCondition* condition )
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  ProposedStep = DBL_MAX;

  if(first) 
  {
    //  checkVolume = histo->CheckVolume();
    //  gasVolume   = histo->GasVolume();
    first = false;
  }

  G4VPhysicalVolume* pv = aTrack.GetVolume();

  if( pv == gasVolume || pv == checkVolume ) ProposedStep = 0.0;

  else if( (aTrack.GetMaterial())->GetDensity() > thDensity && 
           aTrack.GetPosition().z() < 0.0 )                    ProposedStep = MaxChargedStep;

  return ProposedStep;
}

////////////////////////////////////////////////////////////////////////////////////

G4VParticleChange* StepMax::PostStepDoIt(const G4Track& aTrack, const G4Step&)
{
  aParticleChange.Initialize(aTrack);

  if ( ProposedStep == 0.0 ) 
  {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    /*
    if(1 < (Histo::GetPointer())->GetVerbose()) 
    {
      G4cout << "StepMax: " << aTrack.GetDefinition()->GetParticleName()
             << " with energy = " << aTrack.GetKineticEnergy()/MeV
             << " MeV is killed in Check volume at " << aTrack.GetPosition()
             << G4endl;
    }
    */
  }
  return &aParticleChange;
}

////////////////////////////////////////////////////////////////////////////////////////
