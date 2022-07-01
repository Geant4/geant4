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

#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4Region.hh"
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
#include "G4FastSimulationVector.hh"

#include "G4ios.hh"

//-------------------------------------------
//
//        G4FastSimulationManager class
//
//-------------------------------------------

// Class Description:
//  The G4VFastSimulationModel objects are attached to the envelope
// through a G4FastSimulationManager.
//   This object will manage the list of models and will message them
// at tracking time.
//


class G4FastSimulationManager
{
public:  // with description
  //------------------------
  // Constructor/Destructor
  //------------------------
  // Only one Constructor. By default the envelope can
  // be placed n-Times. 
  // If the user is sure that it is placed just one time,
  // the IsUnique flag should be set TRUE to avoid the
  // G4AffineTransform re-calculations each time we reach
  // the envelope.

  G4FastSimulationManager(G4Envelope *anEnvelope,
			  G4bool        IsUnique = FALSE);
  // This is the only constructor. In this constructor you specify 
  // the envelope by giving the G4Region (typedef-ed as G4Envelope)
  // pointer. The G4FastSimulationManager object will bind itself to
  // this envelope and will notify this G4Region to become an envelope.
  // If you know that this region is used for only one logical volume,
  // you can turn the IsUnique boolean to "true" to allow some optimization. 
  //
  // Note that if you choose to use the G4VFastSimulationModel(const G4String&,
  // G4Region*, G4bool) constructor for you model, the G4FastSimulationManager
  // will be constructed using the given G4Region* and G4bool values of the
  // model constructor.
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

public:  // without description  
  // Methods for print/control commands
  void ListTitle() const;
  void ListModels() const;
  void ListModels(const G4ParticleDefinition*) const;
  void ListModels(const G4String& aName) const;
  const G4Envelope* GetEnvelope() const;

  G4VFastSimulationModel* GetFastSimulationModel(const G4String& modelName,
						 const G4VFastSimulationModel* previousFound,
						 bool &foundPrevious) const;

  const std::vector<G4VFastSimulationModel*>& GetFastSimulationModelList() const
  {return ModelList;}

  void FlushModels();

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

  // For management
  G4bool operator == ( const G4FastSimulationManager&) const;

private:
  // Private members :
  G4FastTrack fFastTrack;
  G4FastStep  fFastStep;
  G4VFastSimulationModel* fTriggedFastSimulationModel;
  G4FastSimulationVector <G4VFastSimulationModel> ModelList;
  G4FastSimulationVector <G4VFastSimulationModel> fInactivatedModels;

  G4ParticleDefinition* fLastCrossedParticle;
  G4FastSimulationVector <G4VFastSimulationModel> fApplicableModelList;

  // -- *** depracating, to be dropped @ next major release:
  G4FastSimulationVector <G4Transform3D> GhostPlacements;
};

inline void 
G4FastSimulationManager::AddFastSimulationModel(G4VFastSimulationModel* fsm)
{
  ModelList.push_back(fsm);
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
