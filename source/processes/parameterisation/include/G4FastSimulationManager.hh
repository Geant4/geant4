// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastSimulationManager.hh,v 1.1 1999-01-07 16:14:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//---------------------------------------------------------------
//
//  G4FastSimulationManager.hh
//
//  Description:
//    Manages the Fast Simulation models attached to a envelope.
//
//  History:
//    Oct 97: Verderi && MoraDeFreitas - First Implementation.
//
//---------------------------------------------------------------


#ifndef G4FastSimulationManager_h
#define G4FastSimulationManager_h 1

#include <rw/tpordvec.h>

#include "globals.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VParticleChange.hh"
#include "G4FastTrack.hh"
#include "G4FastStep.hh"
#include "G4VFastSimulationModel.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4ios.hh"

//---------------------------
// For possible future needs:
//---------------------------
typedef G4LogicalVolume G4Envelope;

//-------------------------------------------
//
//        G4FastSimulationManager class
//
//-------------------------------------------
class G4FastSimulationManager
{
public:
  //------------------------
  // Constructor/Destructor
  //------------------------
  // Only one Constructor. By default the envelope can
  // be placed n-Times. 
  // If the user is sure that it'll be placed just one time,
  // the IsUnique flag should be set TRUE to avoid the
  // G4AffineTransform re-calculations each time we reach
  // the envelope.
  G4FastSimulationManager(G4Envelope *anEnvelope,
			  G4bool IsUnique=FALSE);
  ~G4FastSimulationManager();

  // Methods to add/remove models to/from the Model 
  // List.
  //
  void AddFastSimulationModel(G4VFastSimulationModel*);
  void RemoveFastSimulationModel(G4VFastSimulationModel*);

  // Methods to activate/inactivate models from the Model 
  // List.

  G4bool ActivateFastSimulationModel(const G4String&);
  G4bool InActivateFastSimulationModel(const G4String&);

  // Methods to add/remove GhostPlacements to/from the 
  // GhostPlacements List.
  //
  G4Transform3D* AddGhostPlacement(G4RotationMatrix*,
				   const G4ThreeVector&);
  
  G4Transform3D* AddGhostPlacement(G4Transform3D*);
  
  G4bool RemoveGhostPlacement(const G4Transform3D*);
  
  
  // Methods for print/control commands
  void ListTitle() const;
  void ListModels() const;
  void ListModels(const G4ParticleDefinition*) const;
  void ListModels(const G4String& aName) const;
  const G4Envelope* GetEnvelope() const;

  //
  //----------------------------------------------
  // Interface methods for the G4GlobalFastSimulationManager
  //----------------------------------------------
  // Parallel geometry placements
  
  G4bool InsertGhostHereIfNecessary(G4VPhysicalVolume* ,
				    const G4ParticleDefinition&);

  //
  //----------------------------------------------
  // Interface methods for the 
  // G4FastSimulationManagerProcess process.
  //----------------------------------------------
  // Trigger
  G4bool PostStepGetFastSimulationManagerTrigger(const G4Track &,
					 const G4Navigator* a = NULL);
  // DoIt
  void InvokePostStepDoIt();

  // AtRest methods:
  G4bool AtRestGetFastSimulationManagerTrigger(const G4Track &,
					       const G4Navigator* a = NULL);
  void InvokeAtRestDoIt();

  //Final G4VParticleChange* to return to the Stepping.
  G4VParticleChange* GettheParticleChange();

  // For RW management
  G4bool operator == ( const G4FastSimulationManager&) const;

private:
  // Private members :
  G4FastTrack fFastTrack;
  G4FastStep  fFastStep;
  G4VFastSimulationModel* fTriggedFastSimulationModel;
  RWTPtrOrderedVector<G4VFastSimulationModel> ModelList;
  RWTPtrOrderedVector<G4VFastSimulationModel> fInactivatedModels;
  RWTPtrOrderedVector<G4Transform3D> GhostPlacements;

  G4ParticleDefinition* fLastCrossedParticle;
  RWTPtrOrderedVector<G4VFastSimulationModel> fApplicableModelList;
};

inline void 
G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel* fsm)
{
  ModelList.insert(fsm);
  // forces the fApplicableModelList to be rebuild
  fLastCrossedParticle=NULL;
}

inline void 
G4FastSimulationManager::RemoveFastSimulationModel(G4VFastSimulationModel* fsm)
{
  if(!ModelList.remove(fsm)) fInactivatedModels.remove(fsm);
  // forces the fApplicableModelList to be rebuild
  fLastCrossedParticle=NULL;
}

inline G4VParticleChange* G4FastSimulationManager::GettheParticleChange()
{
  return &fFastStep;
}

inline G4bool 
G4FastSimulationManager::operator == (const G4FastSimulationManager& fsm) const
{
  return (this==&fsm) ? true : false;
}

inline const G4Envelope* 
G4FastSimulationManager::GetEnvelope() const
{
  return fFastTrack.GetEnvelope();
}

#endif
