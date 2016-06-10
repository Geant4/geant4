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
// $Id: G4WrapperProcess.cc 80787 2014-05-12 09:06:07Z gcosmo $
//
// 
// ------------------------------------------------------------
//        GEANT 4 class implementation file 
//
// ------------------------------------------------------------
//   New Physics scheme           18 Dec. 1996  H.Kurahige
// ------------------------------------------------------------

#include "G4WrapperProcess.hh"

G4WrapperProcess::G4WrapperProcess(const G4String& aName,
                                         G4ProcessType   aType)
  : G4VProcess(aName,aType), pRegProcess((G4VProcess*)(0))
{
}

G4WrapperProcess::G4WrapperProcess(const G4WrapperProcess& right)
  : G4VProcess(*((G4VProcess*)(&right))), pRegProcess(right.pRegProcess)
{
}

G4WrapperProcess::~G4WrapperProcess()
{
  if (pRegProcess!=0) delete pRegProcess;
}

void G4WrapperProcess::ResetNumberOfInteractionLengthLeft()
{
  pRegProcess->ResetNumberOfInteractionLengthLeft();
}

G4double G4WrapperProcess::
AlongStepGetPhysicalInteractionLength( const G4Track& track,
                                             G4double  previousStepSize,
                                             G4double  currentMinimumStep,
                                             G4double& proposedSafety,
                                             G4GPILSelection* selection     )
{
  return pRegProcess->
         AlongStepGetPhysicalInteractionLength( track,
                                                previousStepSize,
                                                currentMinimumStep,
                                                proposedSafety,
                                                selection     );
}

G4double G4WrapperProcess::
AtRestGetPhysicalInteractionLength( const G4Track& track,
                                          G4ForceCondition* condition )
{
  return pRegProcess->AtRestGetPhysicalInteractionLength( track, condition );
}

G4double G4WrapperProcess::
PostStepGetPhysicalInteractionLength( const G4Track& track,
                                            G4double   previousStepSize,
                                            G4ForceCondition* condition )
{
   return pRegProcess->PostStepGetPhysicalInteractionLength( track,
                                                             previousStepSize,
                                                             condition );
}
      
void G4WrapperProcess::SetProcessManager(const G4ProcessManager* procMan)
{
   pRegProcess->SetProcessManager(procMan); 
}

const G4ProcessManager* G4WrapperProcess::GetProcessManager()
{
  return     pRegProcess->GetProcessManager();
}

G4VParticleChange* G4WrapperProcess::PostStepDoIt( const G4Track& track,
                                                   const G4Step&  stepData )
{
  return     pRegProcess->PostStepDoIt( track, stepData );        
}

G4VParticleChange* G4WrapperProcess::AlongStepDoIt( const G4Track& track,
                                                    const G4Step& stepData )
{
  return     pRegProcess->AlongStepDoIt( track, stepData );        
}
 
G4VParticleChange* G4WrapperProcess::AtRestDoIt( const G4Track& track,
                                                 const G4Step& stepData )
{
  return     pRegProcess->AtRestDoIt( track, stepData );        
}

G4bool G4WrapperProcess::IsApplicable(const G4ParticleDefinition& particle)
{
  return     pRegProcess->IsApplicable(particle);
}

void G4WrapperProcess::BuildPhysicsTable(const G4ParticleDefinition& particle)
{
  return     pRegProcess->BuildPhysicsTable(particle);
}

void G4WrapperProcess::PreparePhysicsTable(const G4ParticleDefinition& particle)
{
  return     pRegProcess->PreparePhysicsTable(particle);
}

G4bool G4WrapperProcess::
StorePhysicsTable(const G4ParticleDefinition* particle,
                  const G4String& directory, 
                        G4bool          ascii)
{
  return pRegProcess->StorePhysicsTable(particle,  directory,  ascii);
} 
 
G4bool G4WrapperProcess::
RetrievePhysicsTable( const G4ParticleDefinition* particle,
                      const G4String& directory, 
                            G4bool          ascii)
{
  return pRegProcess->RetrievePhysicsTable(particle,  directory,  ascii);
}  

void G4WrapperProcess::StartTracking(G4Track* track)
{
  pRegProcess->StartTracking(track);
}

void G4WrapperProcess::EndTracking()
{
  pRegProcess->EndTracking();
}

void   G4WrapperProcess::RegisterProcess(G4VProcess* process)
{
  pRegProcess=process;
  theProcessName += process->GetProcessName();
  theProcessType = process->GetProcessType();
}

const G4VProcess* G4WrapperProcess::GetRegisteredProcess() const
{
  return pRegProcess;
} 

void G4WrapperProcess::SetMasterProcess(G4VProcess* masterP) {
  G4WrapperProcess* master = static_cast<G4WrapperProcess*>(masterP);
  //Cannot use getter because it returns "const" and we do not want
  //to cast away constness explicitly (even if this is the same)
  pRegProcess->SetMasterProcess(master->pRegProcess);
}
