// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastSimulationManager.cc,v 1.3 1999-12-15 14:53:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//$Id:
//---------------------------------------------------------------
//
//  G4FastSimulationManager.cc
//
//  Description:
//    Manages the Fast Simulation models attached to a envelope.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------

#include "G4FastSimulationManager.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4PVPlacement.hh"

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
  // Communicates to the Logical Volume that it becomes a
  // envelope and with this fast simulation manager.
  anEnvelope->BecomeEnvelopeForFastSimulation(this);

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
    fFastTrack.GetEnvelope()->ClearEnvelopeForFastSimulation();
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
  for (iModel=0; iModel<ModelList.length(); iModel++)
    if(ModelList(iModel)->GetName() == aName)
      return true;
  
  // Look for in the fInactivatedModels list, if found insert it back to 
  // the ModelList
  for (iModel=0; iModel<fInactivatedModels.length(); iModel++)
    if(fInactivatedModels(iModel)->GetName() == aName) {
      ModelList.
	insert(fInactivatedModels.removeAt(iModel));
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
  for (G4int iModel=0; iModel<ModelList.length(); iModel++)
    if(ModelList(iModel)->GetName() == aName) {
      fInactivatedModels.
	insert(ModelList.removeAt(iModel));
      // forces the fApplicableModelList to be rebuild
      fLastCrossedParticle=0;
      return true;
    }
  return false;
}

//----------------------------------------
// Methods to add/remove GhostPlacements 
//----------------------------------------

G4Transform3D*
G4FastSimulationManager::AddGhostPlacement(G4RotationMatrix *prot,
					   const G4ThreeVector &tlate)
{
  G4Transform3D* newghostplace;
  if(prot==0) prot = new G4RotationMatrix();
  newghostplace = new G4Transform3D(*prot,tlate);
  AddGhostPlacement(newghostplace);
  return newghostplace;
}

G4Transform3D*
G4FastSimulationManager::AddGhostPlacement(G4Transform3D *trans3d)
{
  GhostPlacements.insert(trans3d);
  G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
    FastSimulationNeedsToBeClosed();  
  return trans3d;
}

G4bool 
G4FastSimulationManager::RemoveGhostPlacement(const G4Transform3D *trans3d)
{
  G4bool found;
  if((found=(GhostPlacements.remove(trans3d) != 0)))
    G4GlobalFastSimulationManager::GetGlobalFastSimulationManager()->
      FastSimulationNeedsToBeClosed();
  return found;
}

//
//-------------------------------------
// Interface trigger method for the 
// G4ParameterisationManagerProcess
//-------------------------------------
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
  G4int iModel;
  
  // If particle type changed re-build the fApplicableModelList.
  if(fLastCrossedParticle!=track.GetDefinition()) {
    fLastCrossedParticle=track.GetDefinition();
    fApplicableModelList.clear();
    // If Model List is empty, do nothing !
    if(ModelList.length()==0) return false;
    for (iModel=0; iModel<ModelList.length(); iModel++)
      if(ModelList(iModel)->IsApplicable(*(track.GetDefinition())))
	fApplicableModelList.insert(ModelList(iModel));
  }

  // If Applicable Model List is empty, do nothing !
  if(fApplicableModelList.length()==0) return false;

  // -- Register current track
  fFastTrack.SetCurrentTrack(track,theNavigator);

  // tests if particle are on the boundary and leaving,
  // in this case do nothing !
  if(fFastTrack.OnTheBoundaryButExiting()) return false;
  
  // Loops on the ModelTrigger() methods
  for (iModel=0; iModel<fApplicableModelList.length(); iModel++)
    
    //---------------------------------------------------
    // Asks the ModelTrigger method if it must be trigged now.
    //---------------------------------------------------
    
    if(fApplicableModelList(iModel)->ModelTrigger(fFastTrack)) {
      //--------------------------------------------------
      // The model will be applied. Initializes the G4FastStep 
      // with the current state of the G4Track and 
      // same usefull parameters.
      // In particular it does SetLocalEnergyDeposit(0.0).
      //--------------------------------------------------	
      fFastStep.Initialize(fFastTrack);
      
      // Keeps the FastSimulationModel pointer to call the
      // DoIt() method.
      fTriggedFastSimulationModel=fApplicableModelList(iModel);
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
  G4int iModel;
  
  // If particle type changed re-build the fApplicableModelList.
  if(fLastCrossedParticle!=track.GetDefinition()) {
    fLastCrossedParticle=track.GetDefinition();
    fApplicableModelList.clear();
    // If Model List is empty, do nothing !
    if(ModelList.length()==0) return false;
    for (iModel=0; iModel<ModelList.length(); iModel++)
      if(ModelList(iModel)->IsApplicable(*(track.GetDefinition())))
	fApplicableModelList.insert(ModelList(iModel));
  }
  
  // If Applicable Model List is empty, do nothing !
  if(fApplicableModelList.length()==0) return false;

  // -- Register current track
  fFastTrack.SetCurrentTrack(track,theNavigator);
  
  // -- (note: compared to the PostStepGetFastSimulationManagerTrigger,
  // --  the test to see if the particle is on the boundary but leaving
  // --  is irrelevant here)
  
  // Loops on the models to see if one of them wants to trigger:
  for (iModel=0; iModel < fApplicableModelList.length(); iModel++)
    if(fApplicableModelList(iModel)->AtRestModelTrigger(fFastTrack))
      {
	fFastStep.Initialize(fFastTrack);
	fTriggedFastSimulationModel=fApplicableModelList(iModel);
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

G4bool 
G4FastSimulationManager::
InsertGhostHereIfNecessary(G4VPhysicalVolume* theClone,
			   const G4ParticleDefinition& theParticle)
{
  G4PVPlacement *GhostPhysical;
  // Not to do if there aren't glost placements
  if(GhostPlacements.length()==0) return false;

  // If there are, verifies if at least one model is applicable
  // for theParticle.
  for (G4int iModel=0; iModel<ModelList.length(); iModel++)
    if(ModelList(iModel)->IsApplicable(theParticle)) {
      // Ok, we find one. Place the ghost(s).
      for (G4int ighost=0; ighost<GhostPlacements.length(); ighost++)
	GhostPhysical=new 
	  G4PVPlacement(*(GhostPlacements(ighost)),
			fFastTrack.GetEnvelope()->GetName(),
			fFastTrack.GetEnvelope(),
			theClone,
			false,0);
      //  And answer true
      return true;
    }
  //otherwise answer false
  return false;  
}

void 
G4FastSimulationManager::ListTitle() const
{
  G4cout << fFastTrack.GetEnvelope()->GetName();
  if(GhostPlacements.length()!=0) G4cout << " (ghost)";
}

void 
G4FastSimulationManager::ListModels() const
{
  G4int iModel;

  G4cout << "Current Models for the ";
  ListTitle();
  G4cout << " Envelope:\n";

  for (iModel=0; iModel<ModelList.length(); iModel++) 
    G4cout << "   " << ModelList(iModel)->GetName() << "\n";

  for (iModel=0; iModel<fInactivatedModels.length(); iModel++)
    G4cout << "   " << fInactivatedModels(iModel)->GetName() 
	   << "(inactivated)\n";
}

void 
G4FastSimulationManager::ListModels(const G4String& aName) const
{
  G4int iModel;
  G4int titled = 0;
  G4ParticleTable* theParticleTable=
    G4ParticleTable::GetParticleTable();
  
  // Active Models
  for (iModel=0; iModel<ModelList.length(); iModel++)
    if(ModelList(iModel)->GetName() == aName ||
       aName == "all" ) {
      if(!(titled++)){
	G4cout << "In the envelope ";
	ListTitle();
	G4cout << ",\n";
      }
      G4cout << "  the model " << ModelList(iModel)->GetName()
	     << " is applicable for :\n     ";
      
      G4int list_started=0;
      for (G4int iParticle=0; iParticle<theParticleTable->entries(); 
	   iParticle++)
	if(ModelList(iModel)->
	   IsApplicable(*(theParticleTable->
			  GetParticle(iParticle)))) {
	  if(list_started++) G4cout << ", ";
	  G4cout << theParticleTable->
	    GetParticle(iParticle)->GetParticleName();
	}
      G4cout <<G4endl;
    }
  
  // Inactive Models
  for (iModel=0; iModel<fInactivatedModels.length(); iModel++)
    if(fInactivatedModels(iModel)->GetName() == aName ||
       aName == "all" ) {
      if(!(titled++)){
	G4cout << "In the envelope ";
	ListTitle();
	G4cout << ",\n";
      }
      G4cout << "  the model " << fInactivatedModels(iModel)->GetName()
	     << " (inactivated) is applicable for :\n     ";
      
      G4int list_started=0;
      for (G4int iParticle=0; iParticle<theParticleTable->entries(); 
	   iParticle++)
	if(fInactivatedModels(iModel)->
	   IsApplicable(*(theParticleTable->
			  GetParticle(iParticle)))) {
	  if(list_started++) G4cout << ", ";
	  G4cout << theParticleTable->
	    GetParticle(iParticle)->GetParticleName();
	}
      G4cout <<G4endl;
    }
}

void 
G4FastSimulationManager::ListModels(const G4ParticleDefinition* aPD) const
{
  G4int iModel;
  G4bool unique=true;
  
  // Active Models
  for (iModel=0; iModel<ModelList.length(); iModel++)
    if(ModelList(iModel)->IsApplicable(*aPD)) {
      G4cout << "Envelope ";
      ListTitle();
      G4cout << ", Model " 
	     << ModelList(iModel)->GetName() 
	     << "." << G4endl;
    }
  // inactive Models
  for (iModel=0; iModel<fInactivatedModels.length(); iModel++)
    if(fInactivatedModels(iModel)->IsApplicable(*aPD)) {
      G4cout << "Envelope ";
      ListTitle();
      G4cout << ", Model " 
	     << fInactivatedModels(iModel)->GetName() 
	     << " (inactivated)." << G4endl;
    }
  
  if(!unique)
    G4cout << "\a\n >>>>>>Warning: two or more Models for the same "
	   << "particle type attached to the same envelope!"
	   << G4endl;
  unique=false;
}
