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
// $Id: G4AdjointProcessEquivalentToDirectProcess.cc 66892 2013-01-17 10:57:59Z gunter $
//
// 
// ------------------------------------------------------------
//        GEANT 4 class implementation file 
//
// Class Description
//
// This class is for adjoint process equivalent to direct process

// ------------------------------------------------------------
//   Created by L.Desorgher          25 Sept. 2009  Inspired from G4WrapperProcess
// ------------------------------------------------------------

#include "G4AdjointProcessEquivalentToDirectProcess.hh"
#include "G4DynamicParticle.hh"
G4AdjointProcessEquivalentToDirectProcess::G4AdjointProcessEquivalentToDirectProcess(const G4String& aName,
										     G4VProcess* aProcess,
										     G4ParticleDefinition* fwd_particle_def)
:G4VProcess(aName)
{  
   theDirectProcess =aProcess;
   theProcessType = theDirectProcess->GetProcessType();
   theFwdParticleDef  = fwd_particle_def;
}


G4AdjointProcessEquivalentToDirectProcess::~G4AdjointProcessEquivalentToDirectProcess()
{
  if (theDirectProcess!=0) delete theDirectProcess;
}

void G4AdjointProcessEquivalentToDirectProcess::ResetNumberOfInteractionLengthLeft()
{
  theDirectProcess->ResetNumberOfInteractionLengthLeft();
}

G4double G4AdjointProcessEquivalentToDirectProcess::
AlongStepGetPhysicalInteractionLength( const G4Track& track,
                                             G4double  previousStepSize,
                                             G4double  currentMinimumStep,
                                             G4double& proposedSafety,
                                             G4GPILSelection* selection     )
{ 
  

  //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts = const_cast<G4DecayProducts*>  (theDynPart->GetPreAssignedDecayProducts());
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  
  //Call the direct process
  //----------------------
  G4double GPIL =  theDirectProcess->
         AlongStepGetPhysicalInteractionLength( track,
                                                previousStepSize,
                                                currentMinimumStep,
                                                proposedSafety,
                                                selection     );
  
  
  //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
  
  return GPIL;
						
}

G4double G4AdjointProcessEquivalentToDirectProcess::
AtRestGetPhysicalInteractionLength( const G4Track& track,
                                          G4ForceCondition* condition )
{ //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts =  const_cast<G4DecayProducts*>  (theDynPart->GetPreAssignedDecayProducts());
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  
  //Call the direct process
  //----------------------

   
  G4double GPIL =  theDirectProcess->AtRestGetPhysicalInteractionLength( track, condition );
  
  //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
  return GPIL;
						
  
}

G4double G4AdjointProcessEquivalentToDirectProcess::
PostStepGetPhysicalInteractionLength( const G4Track& track,
                                            G4double   previousStepSize,
                                            G4ForceCondition* condition )
{
  //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts =  const_cast<G4DecayProducts*>  (theDynPart->GetPreAssignedDecayProducts());
 
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  
  //Call the direct process
  //----------------------

   
  G4double GPIL = theDirectProcess->PostStepGetPhysicalInteractionLength( track,
                                                             previousStepSize,
                                                             condition );
   
  //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
   return GPIL;
  
  				     
}
/*
      
void G4AdjointProcessEquivalentToDirectProcess::SetProcessManager(const G4ProcessManager* procMan)
{
   theDirectProcess->SetProcessManager(procMan); 
}

const G4ProcessManager* G4AdjointProcessEquivalentToDirectProcess::GetProcessManager()
{
  return     theDirectProcess->GetProcessManager();
}
*/
G4VParticleChange* G4AdjointProcessEquivalentToDirectProcess::PostStepDoIt( const G4Track& track,
                                                   const G4Step&  stepData )
{
  //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts =  const_cast<G4DecayProducts*>  (theDynPart->GetPreAssignedDecayProducts());
 
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  
  //Call the direct process
  //----------------------
  
  G4VParticleChange* partChange = theDirectProcess->PostStepDoIt( track, stepData );
  
  
  //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
  return partChange;
  
  
          
}

G4VParticleChange* G4AdjointProcessEquivalentToDirectProcess::AlongStepDoIt( const G4Track& track,
                                                    const G4Step& stepData )
{ 
  //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts =  const_cast<G4DecayProducts*>  (theDynPart->GetPreAssignedDecayProducts());
 
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  
  //Call the direct process
  //----------------------
  G4VParticleChange* partChange =theDirectProcess->AlongStepDoIt( track, stepData );
  
  //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
  return partChange;        
}
 
G4VParticleChange* G4AdjointProcessEquivalentToDirectProcess::AtRestDoIt( const G4Track& track,
                                                 const G4Step& stepData )
{
  //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track.GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts =  const_cast<G4DecayProducts*>  (theDynPart->GetPreAssignedDecayProducts());
 
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  
  //Call the direct process
  //----------------------
  G4VParticleChange* partChange =theDirectProcess->AtRestDoIt( track, stepData );
  
  //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
   return partChange; 
  
       
}

G4bool G4AdjointProcessEquivalentToDirectProcess::IsApplicable(const G4ParticleDefinition&)
{
  return     theDirectProcess->IsApplicable(*theFwdParticleDef);
}

void G4AdjointProcessEquivalentToDirectProcess::BuildPhysicsTable(const G4ParticleDefinition& )
{
  return     theDirectProcess->BuildPhysicsTable(*theFwdParticleDef);
}

void G4AdjointProcessEquivalentToDirectProcess::PreparePhysicsTable(const G4ParticleDefinition& )
{
  return     theDirectProcess->PreparePhysicsTable(*theFwdParticleDef);
}

G4bool G4AdjointProcessEquivalentToDirectProcess::
StorePhysicsTable(const G4ParticleDefinition* ,
                  const G4String& directory, 
                        G4bool          ascii)
{
  return theDirectProcess->StorePhysicsTable(theFwdParticleDef,  directory,  ascii);
} 
 
G4bool G4AdjointProcessEquivalentToDirectProcess::
RetrievePhysicsTable( const G4ParticleDefinition* ,
                      const G4String& directory, 
                            G4bool          ascii)
{
  return theDirectProcess->RetrievePhysicsTable(theFwdParticleDef,  directory,  ascii);
}  

void G4AdjointProcessEquivalentToDirectProcess::StartTracking(G4Track* track)
{
  //Change the particle definition to the direct one
  //------------------------------------------------
  G4DynamicParticle* theDynPart = const_cast<G4DynamicParticle*> (track->GetDynamicParticle());
  G4ParticleDefinition* adjPartDef = theDynPart->GetDefinition();
  
  G4DecayProducts* decayProducts =  const_cast<G4DecayProducts*> (theDynPart->GetPreAssignedDecayProducts());
  theDynPart->SetPreAssignedDecayProducts((G4DecayProducts*)(0));
  theDynPart->SetDefinition(theFwdParticleDef);
  
  theDirectProcess->StartTracking(track);
  
   //Restore the adjoint particle definition to the direct one
  //------------------------------------------------
  theDynPart->SetDefinition(adjPartDef);
  theDynPart->SetPreAssignedDecayProducts(decayProducts);
  
  
  return; 
  
}

void G4AdjointProcessEquivalentToDirectProcess::EndTracking()
{
  theDirectProcess->EndTracking(); 
}

