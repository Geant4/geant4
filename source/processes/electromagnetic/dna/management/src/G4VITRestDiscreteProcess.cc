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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4VITRestDiscreteProcess.hh"

G4VITRestDiscreteProcess::G4VITRestDiscreteProcess(const G4String& aName,
                                                   G4ProcessType aType)
    : G4VITProcess(aName, aType)
{
    enableAlongStepDoIt = false;
}

G4VITRestDiscreteProcess::~G4VITRestDiscreteProcess() = default;

G4double G4VITRestDiscreteProcess::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                                        G4double previousStepSize,
                                                                        G4ForceCondition* condition)
{
    if ((previousStepSize < 0.0) || (fpState->theNumberOfInteractionLengthLeft
                                     <= 0.0))
    {
        // beggining of tracking (or just after DoIt of this process)
        ResetNumberOfInteractionLengthLeft();
    }
    else if (previousStepSize > 0.0)
    {
        // subtract NumberOfInteractionLengthLeft
        SubtractNumberOfInteractionLengthLeft(previousStepSize);
    }
    else
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
    if (fpState->currentInteractionLength < DBL_MAX)
    {
        value = fpState->theNumberOfInteractionLengthLeft * (fpState->currentInteractionLength);
    }
    else
    {
        value = DBL_MAX;
    }
#ifdef G4VERBOSE
    if (verboseLevel > 1)
    {
        G4cout << "G4VITRestDiscreteProcess::PostStepGetPhysicalInteractionLength ";
        G4cout << "[ " << GetProcessName() << "]" << G4endl;
        track.GetDynamicParticle()->DumpInfo();
        G4cout << " in Material  " << track.GetMaterial()->GetName() << G4endl;
        G4cout << "InteractionLength= " << value / CLHEP::cm << "[cm] " << G4endl;
    }
#endif
    return value;
}

G4VParticleChange* G4VITRestDiscreteProcess::PostStepDoIt(const G4Track&,
                                                          const G4Step&)
{
//  reset NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}

G4double G4VITRestDiscreteProcess::AtRestGetPhysicalInteractionLength(const G4Track& track,
                                                                      G4ForceCondition* condition)
{
    // beggining of tracking
    ResetNumberOfInteractionLengthLeft();

    // condition is set to "Not Forced"
    *condition = NotForced;

    // get mean life time
    fpState->currentInteractionLength = GetMeanLifeTime(track, condition);

#ifdef G4VERBOSE
    if ((fpState->currentInteractionLength < 0.0) || (verboseLevel > 2))
    {
        G4cout << "G4VITRestDiscreteProcess::AtRestGetPhysicalInteractionLength ";
        G4cout << "[ " << GetProcessName() << "]" << G4endl;
        track.GetDynamicParticle()->DumpInfo();
        G4cout << " in Material  " << track.GetMaterial()->GetName() << G4endl;
        G4cout << "MeanLifeTime = " << fpState->currentInteractionLength / CLHEP::ns
               << "[ns]" << G4endl;
    }
#endif

    return fpState->theNumberOfInteractionLengthLeft
           * (fpState->currentInteractionLength);
}

G4VParticleChange* G4VITRestDiscreteProcess::AtRestDoIt(const G4Track&,
                                                        const G4Step&)
{
//  clear NumberOfInteractionLengthLeft
    ClearNumberOfInteractionLengthLeft();

    return pParticleChange;
}

