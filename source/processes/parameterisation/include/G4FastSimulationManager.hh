// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FastSimulationManager.hh,v 1.5 1999-12-15 14:53:45 gunter Exp $
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

#include "g4rw/tpordvec.h"

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

// Class Description:
//  The G4VFastSimulationModel objects are attached to the envelope through a G4FastSimulationManager.
//   This object will manage the list of models and will message them at tracking time.
//


class G4FastSimulationManager
{
public:  // with description
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
  // This is the only constructor. In this constructor you specify the envelope by giving the 
  // G4LogicalVolume pointer. The G4FastSimulationManager object will bind itself to this envelope 
  // and will notify this G4LogicalVolume to become an envelope. If you know that this volume is
  // placed only once, you can turn the IsUnique boolean to "true" to allow some optimization. 
  //
  // Note that if you choose to use the G4VFastSimulationModel(const G4String&, G4LogicalVolume*, 
  // G4bool) constructor for you model, the G4FastSimulationManager will be constructed using the 
  // given G4LogicalVolume* and G4bool values of the model constructor.
  //

public:  // without description
  ~G4FastSimulationManager();


public:  // with description
  // Methods to add/remove models to/from the Model 
  // List.
  //
  void AddFastSimulationModel(G4VFastSimulationModel*);
  // Add a model to the Model List.

  void RemoveFastSimulationModel(G4VFastSimulationModel*);
  // Remove a model from the Model List.

  // Methods to activate/inactivate models from the Model 
  // List.

  G4bool ActivateFastSimulationModel(const G4String&);
  // Activate a model in the Model List.

  G4bool InActivateFastSimulationModel(const G4String&);
  // Inactivate a model in the Model List.

  // Methods to add/remove GhostPlacements to/from the 
  // GhostPlacements List.
  //
  G4Transform3D* AddGhostPlacement(G4RotationMatrix*,
				   const G4ThreeVector&);
  // Flag that the envelope is a ghost volume giving its global placement, where the rotation matrix
  // and the translatation vector of 3D transformation describe the placement relative to the world 
  // coordinates. 

  G4Transform3D* AddGhostPlacement(G4Transform3D*);
  // The same but using a G4Transform3D.
  
  G4bool RemoveGhostPlacement(const G4Transform3D*);
  // Removes a Ghost placement.

public:  // without description  
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
					 const G4Navigator* a = 0);
  // DoIt
  G4VParticleChange* InvokePostStepDoIt();

  // AtRest methods:
  G4bool AtRestGetFastSimulationManagerTrigger(const G4Track &,
					       const G4Navigator* a = 0);
  G4VParticleChange*  InvokeAtRestDoIt();

  // For RW management
  G4bool operator == ( const G4FastSimulationManager&) const;

private:
  // Private members :
  G4FastTrack fFastTrack;
  G4FastStep  fFastStep;
  G4VFastSimulationModel* fTriggedFastSimulationModel;
  G4RWTPtrOrderedVector<G4VFastSimulationModel> ModelList;
  G4RWTPtrOrderedVector<G4VFastSimulationModel> fInactivatedModels;
  G4RWTPtrOrderedVector<G4Transform3D> GhostPlacements;

  G4ParticleDefinition* fLastCrossedParticle;
  G4RWTPtrOrderedVector<G4VFastSimulationModel> fApplicableModelList;
};

inline void 
G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel* fsm)
{
  ModelList.insert(fsm);
  // forces the fApplicableModelList to be rebuild
  fLastCrossedParticle = 0;
}

inline void 
G4FastSimulationManager::RemoveFastSimulationModel(G4VFastSimulationModel* fsm)
{
  if(!ModelList.remove(fsm)) fInactivatedModels.remove(fsm);
  // forces the fApplicableModelList to be rebuild
  fLastCrossedParticle = 0;
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
