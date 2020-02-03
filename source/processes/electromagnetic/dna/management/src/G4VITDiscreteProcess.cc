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
//
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
//	History: first implementation, based on object model of
//	2nd December 1995, G.Cosmo
// --------------------------------------------------------------
//   New Physics scheme           8 Jan. 1997  H.Kurahige
// ------------------------------------------------------------

#include "G4VITDiscreteProcess.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

G4VITDiscreteProcess::G4VITDiscreteProcess() :
    G4VITProcess("No Name Discrete Process")
{
  G4Exception("G4VDiscreteProcess::G4VDiscreteProcess()",
              "ProcMan102",
              JustWarning,
              "Default constructor is called");
}

//------------------------------------------------------------------------------

G4VITDiscreteProcess::G4VITDiscreteProcess(const G4String& aName,
                                           G4ProcessType aType) :
    G4VITProcess(aName, aType)
{
  enableAtRestDoIt = false;
  enableAlongStepDoIt = false;
}

//------------------------------------------------------------------------------

G4VITDiscreteProcess::~G4VITDiscreteProcess()
{
}

//------------------------------------------------------------------------------

G4VITDiscreteProcess::G4VITDiscreteProcess(G4VITDiscreteProcess& right) :
    G4VITProcess(right)
{
}

//------------------------------------------------------------------------------

G4double
G4VITDiscreteProcess::
PostStepGetPhysicalInteractionLength(const G4Track& track,
                                     G4double previousStepSize,
                                     G4ForceCondition* condition)
{
  if((previousStepSize < 0.0)
      || (fpState->theNumberOfInteractionLengthLeft <= 0.0))
  {
    // beginning of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else if(previousStepSize > 0.0)
  {
    // subtract NumberOfInteractionLengthLeft
    SubtractNumberOfInteractionLengthLeft(previousStepSize);
  } else
  {
    // zero step
    //  DO NOTHING
  }

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean free path
  fpState->currentInteractionLength = GetMeanFreePath(track,
                                             previousStepSize,
                                             condition);

  G4double value;
  if(fpState->currentInteractionLength < DBL_MAX)
  {
    value = fpState->theNumberOfInteractionLengthLeft
            * fpState->currentInteractionLength;
  } else
  {
    value = DBL_MAX;
  }
#ifdef G4VERBOSE
  if(verboseLevel > 1)
  {
    G4cout << "G4VDiscreteProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" << G4endl;
    track.GetDynamicParticle()->DumpInfo();
    G4cout << " in Material  " << track.GetMaterial()->GetName() << G4endl;
    G4cout << "InteractionLength= " << value / cm << "[cm] " << G4endl;
  }
#endif
  return value;
}

//------------------------------------------------------------------------------

G4VParticleChange* G4VITDiscreteProcess::PostStepDoIt(const G4Track&,
                                                      const G4Step&)
{
//  clear NumberOfInteractionLengthLeft
  ClearNumberOfInteractionLengthLeft();

  return pParticleChange;
}
