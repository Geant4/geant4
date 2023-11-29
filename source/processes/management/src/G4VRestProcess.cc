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
// G4VRestProcess class implementation
//
// Authors:
// - 2 December 1995, G.Cosmo - First implementation, based on object model
// - 18 December 1996, H.Kurashige - New Physics scheme
// --------------------------------------------------------------------

#include "G4VRestProcess.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4MaterialTable.hh"
#include "G4VParticleChange.hh"

// --------------------------------------------------------------------
G4VRestProcess::G4VRestProcess()
  : G4VProcess("No Name Rest Process") 
{
  G4Exception("G4VRestProcess::G4VRestProcess()", "ProcMan102",
              JustWarning, "Default constructor is called");
}

// --------------------------------------------------------------------
G4VRestProcess::G4VRestProcess(const G4String& aName, G4ProcessType aType)
  : G4VProcess(aName, aType)
{
  enableAlongStepDoIt = false;
  enablePostStepDoIt = false;
}

// --------------------------------------------------------------------
G4VRestProcess::~G4VRestProcess()
{
}

// --------------------------------------------------------------------
G4VRestProcess::G4VRestProcess(G4VRestProcess& right)
  : G4VProcess(right)
{
}

// --------------------------------------------------------------------
G4double G4VRestProcess::
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
    G4cout << "G4VRestProcess::AtRestGetPhysicalInteractionLength() - ";
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
G4VParticleChange* G4VRestProcess::AtRestDoIt( const G4Track&,
                                               const G4Step& )
{
  ClearNumberOfInteractionLengthLeft();

  return pParticleChange;
}
