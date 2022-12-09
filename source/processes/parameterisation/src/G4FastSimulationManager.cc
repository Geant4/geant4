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
//---------------------------------------------------------------
//
//  G4FastSimulationManager.cc
//
//  Description:
//    Manages the Fast Simulation models attached to a envelope.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//    ...
//    May 07: Move to parallel world scheme
//
//---------------------------------------------------------------

#include "G4FastSimulationManager.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"

// --------------------------------------------------
// Constructor with envelope and IsUnique flag :
// --------------------------------------------------
//
G4FastSimulationManager::
G4FastSimulationManager(G4Envelope *anEnvelope,
			G4bool IsUnique) :
  fFastTrack(anEnvelope,IsUnique),fTriggedFastSimulationModel(0),
  fLastCrossedParticle(0)
{
  // Communicates to the region that it becomes a
  // envelope and with this fast simulation manager.
  anEnvelope->SetFastSimulationManager(this);

  // Add itself to the GlobalFastSimulationManager 
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
    AddFastSimulationManager(this);
}

// -----------
// Destructor:
// -----------
G4FastSimulationManager::~G4FastSimulationManager()
{
  //
  // Check out the Envelope about this pointer. If in use, 
  // resets the Logical Volume IsEnvelope flag to avoid clash.
  //
  if(fFastTrack.GetEnvelope()->GetFastSimulationManager()==this)
    fFastTrack.GetEnvelope()->ClearFastSimulationManager();
  // Remove itself from the GlobalFastSimulationManager 
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
    RemoveFastSimulationManager(this);
}

// ---------------------------------------
// Methods to activate/inactivate models
//----------------------------------------

G4bool
G4FastSimulationManager::ActivateFastSimulationModel(const G4String& aName) 
{
  G4int iModel;

  // If the model is already active, do nothing.
  for (iModel=0; iModel<(G4int)ModelList.size(); ++iModel)
    if(ModelList[iModel]->GetName() == aName)
      return true;
  
  // Look for in the fInactivatedModels list, if found push_back it back to 
  // the ModelList
  for (iModel=0; iModel<(G4int)fInactivatedModels.size(); ++iModel)
    if(fInactivatedModels[iModel]->GetName() == aName) {
      ModelList.
	push_back (fInactivatedModels.removeAt(iModel));
      // forces the fApplicableModelList to be rebuild
      fLastCrossedParticle=0;
      return true;
    }
  return false;
}

G4bool
G4FastSimulationManager::InActivateFastSimulationModel(const G4String& aName)
{
  // Look for in the ModelList, if found remove from it and keep the pointer 
  // on the fInactivatedModels list.
  for (G4int iModel=0; iModel<(G4int)ModelList.size(); ++iModel)
    if(ModelList[iModel]->GetName() == aName) {
      fInactivatedModels.push_back (ModelList.removeAt(iModel));
      // forces the fApplicableModelList to be rebuild
      fLastCrossedParticle=0;
      return true;
    }
  return false;
}

G4VFastSimulationModel* 
G4FastSimulationManager::GetFastSimulationModel(const G4String& modelName,
						const G4VFastSimulationModel* previousFound,
						bool &foundPrevious) const
{
  G4VFastSimulationModel* model = 0;
  for (std::size_t iModel=0; iModel<ModelList.size(); ++iModel)
    {
      if(ModelList[iModel]->GetName() == modelName)
	{
	  if (previousFound == 0)
	    {
	      model = ModelList[iModel];
	      break;
	    }
	  else
	    {
	      if (ModelList[iModel] == previousFound)
		{
		  foundPrevious = true;
		  continue;
		}
	      if (foundPrevious)
		{
		  model = ModelList[iModel];
		  break;
		}
	    }
	}
    }
  return model;
}

void G4FastSimulationManager::FlushModels()
{
  for (std::size_t iModel=0; iModel<ModelList.size(); ++iModel)
  {
    ModelList[iModel]->Flush();
  }  
}



//------------------------------------------------------------------
// Interface trigger method for the G4ParameterisationManagerProcess
//------------------------------------------------------------------
//   G4bool GetFastSimulationManagerTrigger(const G4Track &);
//
//    This method is used to interface the G4FastSimulationManagerProcess
//    with the user Fast Simulation Models. It's called when the particle 
//    is inside the envelope.
//
//    It :
//
//      1) initialises the private members (fFastTrack and so
//         on);
//      2) loops on the IsApplicable() methods to find out the
//         ones should be applied.
//      2) for these, loops on the ModelTrigger() methods to find out 
//         perhaps one that must be applied just now.
//
//    If the a Fast Simulation Model is triggered then it returns 
//    true, false otherwise.
//
//-----------------------------------------------------------
G4bool 
G4FastSimulationManager::
PostStepGetFastSimulationManagerTrigger(const G4Track& track,
					const G4Navigator* theNavigator)
{
  std::size_t iModel;
  
  // If particle type changed re-build the fApplicableModelList.
  if(fLastCrossedParticle!=track.GetDefinition()) {
    fLastCrossedParticle=track.GetDefinition();
    fApplicableModelList.clear();
    // If Model List is empty, do nothing !
    if(ModelList.size()==0) return false;
    for (iModel=0; iModel<ModelList.size(); ++iModel)
      if(ModelList[iModel]->IsApplicable(*(track.GetDefinition())))
	fApplicableModelList.push_back (ModelList[iModel]);
  }

  // If Applicable Model List is empty, do nothing !
  if(fApplicableModelList.size()==0) return false;

  // -- Register current track
  fFastTrack.SetCurrentTrack(track,theNavigator);

  // tests if particle are on the boundary and leaving,
  // in this case do nothing !
  if(fFastTrack.OnTheBoundaryButExiting()) return false;
  
  // Loops on the ModelTrigger() methods
  for (iModel=0; iModel<fApplicableModelList.size(); ++iModel)
    
    //---------------------------------------------------
    // Asks the ModelTrigger method if it must be trigged now.
    //---------------------------------------------------
    
    if(fApplicableModelList[iModel]->ModelTrigger(fFastTrack)) {
      //--------------------------------------------------
      // The model will be applied. Initializes the G4FastStep 
      // with the current state of the G4Track and 
      // same usefull parameters.
      // In particular it does SetLocalEnergyDeposit(0.0).
      //--------------------------------------------------	
      fFastStep.Initialize(fFastTrack);
      
      // Keeps the FastSimulationModel pointer to call the
      // DoIt() method.
      fTriggedFastSimulationModel=fApplicableModelList[iModel];
      return true;
    }

  //--------------------------------------------
  // Nobody asks to gain control, returns false
  //--------------------------------------------
  return false;
}

G4VParticleChange* G4FastSimulationManager::InvokePostStepDoIt() 
{
  //  const G4FastTrack& parFastTrack=fFastTrack;
  fTriggedFastSimulationModel->DoIt(fFastTrack,fFastStep);
  return &fFastStep;
}

// -------------------------------------------------------------
// -- Mostly the same as above, in the case of AtRest particles:
// -------------------------------------------------------------
G4bool 
G4FastSimulationManager::AtRestGetFastSimulationManagerTrigger(const G4Track& track,
							       const G4Navigator* theNavigator)
{
  std::size_t iModel;
  
  // If particle type changed re-build the fApplicableModelList.
  if(fLastCrossedParticle!=track.GetDefinition()) {
    fLastCrossedParticle=track.GetDefinition();
    fApplicableModelList.clear();
    // If Model List is empty, do nothing !
    if(ModelList.size()==0) return false;
    for (iModel=0; iModel<ModelList.size(); ++iModel)
      if(ModelList[iModel]->IsApplicable(*(track.GetDefinition())))
	fApplicableModelList.push_back (ModelList[iModel]);
  }
  
  // If Applicable Model List is empty, do nothing !
  if(fApplicableModelList.size()==0) return false;

  // -- Register current track
  fFastTrack.SetCurrentTrack(track,theNavigator);
  
  // -- (note: compared to the PostStepGetFastSimulationManagerTrigger,
  // --  the test to see if the particle is on the boundary but leaving
  // --  is irrelevant here)
  
  // Loops on the models to see if one of them wants to trigger:
  for (iModel=0; iModel < fApplicableModelList.size(); ++iModel)
    if(fApplicableModelList[iModel]->AtRestModelTrigger(fFastTrack))
      {
	fFastStep.Initialize(fFastTrack);
	fTriggedFastSimulationModel=fApplicableModelList[iModel];
	return true;
      }
  
  //--------------------------------------------
  // Nobody asks to gain control, returns false
  //--------------------------------------------
  return false;
}

G4VParticleChange* G4FastSimulationManager::InvokeAtRestDoIt() 
{
  fTriggedFastSimulationModel->AtRestDoIt(fFastTrack,fFastStep);
  return &fFastStep;
}

void 
G4FastSimulationManager::ListTitle() const
{
  G4cout << fFastTrack.GetEnvelope()->GetName();
  //  if(GhostPlacements.size()!=0) G4cout << " (ghost)";
  if (fFastTrack.GetEnvelope()->GetWorldPhysical() == G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()) G4cout << " (mass geom.)";
  else G4cout << " (// geom.)";
																			  
}

void 
G4FastSimulationManager::ListModels() const
{
  std::size_t iModel;

  G4cout << "Current Models for the ";
  ListTitle();
  G4cout << " envelope:\n";

  for (iModel=0; iModel<ModelList.size(); ++iModel) 
    G4cout << "   " << ModelList[iModel]->GetName() << "\n";

  for (iModel=0; iModel<fInactivatedModels.size(); ++iModel)
    G4cout << "   " << fInactivatedModels[iModel]->GetName() 
	   << "(inactivated)\n";
}

void G4FastSimulationManager::ListModels(const G4String& modelName) const
{
  std::size_t iModel;
  G4int titled = 0;
  G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  
  // Active Models
  for ( iModel=0; iModel<ModelList.size(); ++iModel )
    if( ModelList[iModel]->GetName() == modelName || modelName == "all" )
      {
	if( !(titled++) )
	  {
	    G4cout << "In the envelope ";
	    ListTitle();
	    G4cout << ",\n";
	  }
	G4cout << "  the model " << ModelList[iModel]->GetName()
	       << " is applicable for :\n     ";
	
	G4int list_started=0;
	for ( G4int iParticle = 0; iParticle<theParticleTable->entries(); iParticle++)
	  if( ModelList[iModel] -> IsApplicable( *(theParticleTable->GetParticle(iParticle))) )
	    {
	      if(list_started++) G4cout << ", ";
	      G4cout << theParticleTable->
		GetParticle(iParticle)->GetParticleName();
	    }
	G4cout <<G4endl;
      }
  
  // Inactive Models
  for (iModel=0; iModel<fInactivatedModels.size(); ++iModel)
    if(fInactivatedModels[iModel]->GetName() == modelName || modelName == "all" )
      {
	if( !(titled++) )
	  {
	    G4cout << "In the envelope ";
	    ListTitle();
	    G4cout << ",\n";
	  }
	G4cout << "  the model " << fInactivatedModels[iModel]->GetName()
	       << " (inactivated) is applicable for :\n     ";
	
	G4int list_started=0;
	for ( G4int iParticle=0; iParticle<theParticleTable->entries(); iParticle++ )
	  if( fInactivatedModels[iModel] -> IsApplicable( *(theParticleTable->GetParticle(iParticle))) )
	    {
	      if(list_started++) G4cout << ", ";
	      G4cout << theParticleTable->
		GetParticle(iParticle)->GetParticleName();
	    }
	G4cout <<G4endl;
      }
}

void G4FastSimulationManager::ListModels(const G4ParticleDefinition* particleDefinition) const
{
  std::size_t iModel;
  G4bool unique = true;
  
  // Active Models
  for ( iModel=0; iModel<ModelList.size(); ++iModel )
    if ( ModelList[iModel]->IsApplicable(*particleDefinition) )
      {
	G4cout << "Envelope ";
	ListTitle();
	G4cout << ", Model " 
	       << ModelList[iModel]->GetName() 
	       << "." << G4endl;
	// -- Verify unicity of model attached to particleDefinition:
	for ( auto jModel = iModel + 1; jModel < ModelList.size(); jModel++ )
	  if ( ModelList[jModel]->IsApplicable(*particleDefinition) ) unique = false;
      }
  
  // Inactive Models
  for ( iModel=0; iModel<fInactivatedModels.size(); ++iModel )
    if( fInactivatedModels[iModel]->IsApplicable(*particleDefinition) )
      {
	G4cout << "Envelope ";
	ListTitle();
	G4cout << ", Model " 
	       << fInactivatedModels[iModel]->GetName() 
	       << " (inactivated)." << G4endl;
      }
  
  if( !unique )
    {
      G4ExceptionDescription ed;
      ed << "Two or more active Models are available for the same particle type, in the same envelope/region." << G4endl;
      G4Exception("G4FastSimulationManager::ListModels(const G4ParticleDefinition* particleDefinition) const",
		  "FastSim001",
		  JustWarning, ed,
		  "Models risk to exclude each other.");
    }
  unique=false;
}
