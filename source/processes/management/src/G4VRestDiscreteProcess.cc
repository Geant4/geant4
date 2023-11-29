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
// G4VRestDiscreteProcess class implementation
//
// Authors:
// - 2 December 1995, G.Cosmo - First implementation, based on object model
// - 8 January 1997, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------

#include "G4VRestDiscreteProcess.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"              
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

// --------------------------------------------------------------------
G4VRestDiscreteProcess::G4VRestDiscreteProcess()
  : G4VProcess("No Name Discrete Process") 
{
  G4Exception("G4VRestDiscreteProcess::G4VRestDiscreteProcess", "ProcMan102",
              JustWarning, "Default constructor is called");
}

// --------------------------------------------------------------------
G4VRestDiscreteProcess::
G4VRestDiscreteProcess(const G4String& aName, G4ProcessType aType)
  : G4VProcess(aName, aType)
{
  enableAlongStepDoIt = false;
}

// --------------------------------------------------------------------
G4VRestDiscreteProcess::~G4VRestDiscreteProcess()
{
}

// --------------------------------------------------------------------
G4VRestDiscreteProcess::
G4VRestDiscreteProcess(G4VRestDiscreteProcess& right)
  : G4VProcess(right)
{
}

// --------------------------------------------------------------------
G4double G4VRestDiscreteProcess::
PostStepGetPhysicalInteractionLength( const G4Track& track,
                                      G4double   previousStepSize,
                                      G4ForceCondition* condition )
{
  if ( (previousStepSize < 0.0) || (theNumberOfInteractionLengthLeft<=0.0))
  {
    // beginning of tracking (or just after DoIt() of this process)
    ResetNumberOfInteractionLengthLeft();
  }
  else if ( previousStepSize > 0.0)
  {
    // subtract NumberOfInteractionLengthLeft 
    SubtractNumberOfInteractionLengthLeft(previousStepSize);
  }
  else
  {
    // zero step
    // DO NOTHING
  }

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean free path
  currentInteractionLength = GetMeanFreePath(track,previousStepSize,condition);

  G4double value;
  if (currentInteractionLength < DBL_MAX)
  {
    value = theNumberOfInteractionLengthLeft * currentInteractionLength;
  }
  else
  {
    value = DBL_MAX;
  }
#ifdef G4VERBOSE
  if (verboseLevel>1)
  {
    G4double length = (value < DBL_MAX) ? value/cm : value;
    G4cout << "G4VRestDiscreteProcess::PostStepGetPhysicalInteractionLength() - ";
    G4cout << "[ " << GetProcessName() << "]" << G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " <<  track.GetMaterial()->GetName() << G4endl;
    G4cout << "InteractionLength= " << length <<"[cm] " << G4endl;
  }
#endif
  return value;
}

// --------------------------------------------------------------------
G4VParticleChange*
G4VRestDiscreteProcess::PostStepDoIt( const G4Track& ,
                                      const G4Step& )
{ 
  ClearNumberOfInteractionLengthLeft();

  return pParticleChange;
}

// --------------------------------------------------------------------
G4double G4VRestDiscreteProcess::
AtRestGetPhysicalInteractionLength( const G4Track& track,
                                    G4ForceCondition* condition )
{
  // beginning of tracking 
  ResetNumberOfInteractionLengthLeft();

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean life time
  currentInteractionLength = GetMeanLifeTime(track, condition);

  G4double time = (currentInteractionLength < DBL_MAX) ?
    theNumberOfInteractionLengthLeft * currentInteractionLength : DBL_MAX;
 
#ifdef G4VERBOSE
  if ((currentInteractionLength <0.0) || (verboseLevel>2))
  {
    G4double t = (currentInteractionLength < DBL_MAX) ?
      currentInteractionLength/ns : DBL_MAX; 
    G4cout << "G4VRestDiscreteProcess::AtRestGetPhysicalInteractionLength() - ";
    G4cout << "[ " << GetProcessName() << "]" << G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() << G4endl;
    G4cout << "MeanLifeTime = " << t << " [ns]"
           << G4endl;
  }
#endif

  return time;
}

// --------------------------------------------------------------------
G4VParticleChange*
G4VRestDiscreteProcess::AtRestDoIt( const G4Track&,
                                    const G4Step& )
{
  ClearNumberOfInteractionLengthLeft();

  return pParticleChange;
}
