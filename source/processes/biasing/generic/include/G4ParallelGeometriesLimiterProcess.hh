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
// $Id: G4ParallelGeometriesLimiterProcess.hh $
//
// 
//---------------------------------------------------------------
//
//  G4ParallelGeometriesLimiterProcess.hh
//
//  Description:
//    Process dedicated to limiting the step on boudaries of
//    parallel geometries used for the generic biasing.
//
//  History:
//    Sep 16: Creation.
//
//---------------------------------------------------------------


#ifndef G4ParallelGeometriesLimiterProcess_hh
#define G4ParallelGeometriesLimiterProcess_hh

#include "globals.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChangeForNothing.hh"
#include "G4FieldTrack.hh"
class G4PathFinder;
class G4TransportationManager;


// Class Description:
// -- G4VProcess limiting the step on parallel geometries used by generic biasing.


class G4ParallelGeometriesLimiterProcess : public G4VProcess
{
public:
  // --------------------------
  // -- Constructor/Destructor:
  // --------------------------
  G4ParallelGeometriesLimiterProcess(const G4String& processName  = "biasLimiter");

public:
  virtual ~G4ParallelGeometriesLimiterProcess()
  {}
  
  // ----------------------------------------------------
  // -- Registration / deregistration of parallel worlds.
  // -- Access to the list of registered parallel worlds.
  // ----------------------------------------------------
  void    AddParallelWorld(const G4String& parallelWorldName);
  void RemoveParallelWorld(const G4String& parallelWorldName);
  // -- The list of registered parallel worlds:
  const std::vector< G4VPhysicalVolume* >&       GetParallelWorlds() const { return fParallelWorlds; }

  // -- Get the parallel world volume index in the list of world volumes handled
  // -- by the process. This index can then be used to access current volume and step
  // -- limitation status below.
  // -- If the world passed is unknown to the process, -1 is returned.
  G4int GetParallelWorldIndex( const G4VPhysicalVolume* parallelWorld     ) const;
  G4int GetParallelWorldIndex( G4String                 parallelWorldName ) const;

  
  // ---------------------
  // -- Active navigators:
  // ---------------------
  // -- The list of navigators handled by the process:
  const std::vector< G4Navigator* >&   GetActiveNavigators()                   const { return fParallelWorldNavigators;                     }
  // -- The navigator used for the passed parallel world index (obtained with GetParallelWorldIndex(...) above)
  // -- Note that no boundary checks are done on the index passed.
  const G4Navigator*                          GetNavigator( G4int worldIndex ) const { return fParallelWorldNavigators[size_t(worldIndex)]; }
  // ---------------------------------------------------
  // -- Previous and current steps geometry information:
  // ---------------------------------------------------
  // -- The "switch" between the previous and current step is done in the PostStepGPIL
  // -- The update on the current step is done:
  // --    - in the PostStepGPIL for the volumes
  // --    - in the AlongStepGPIL for the step limitations
  // --
  // -- The list of previous step and current step volumes:
  const std::vector< const G4VPhysicalVolume* >&  GetCurrentVolumes()                   const { return  fCurrentVolumes;                              }
  const std::vector< const G4VPhysicalVolume* >& GetPreviousVolumes()                   const { return fPreviousVolumes;                              }
  // -- The current and previous volume for the passed parallel world index (obtained with GetParallelWorldIndex(...) above)
  // -- Note that no boundary checks are done on the index passed.
  const G4VPhysicalVolume*                         GetCurrentVolume( G4int worldIndex ) const { return  fCurrentVolumes[size_t(worldIndex)];          }
  const G4VPhysicalVolume*                        GetPreviousVolume( G4int worldIndex ) const { return fPreviousVolumes[size_t(worldIndex)];          }
  // -- Flags telling about step limitation in previous and current step:
  const std::vector< G4bool >&                        GetIsLimiting()                   const { return  fParallelWorldIsLimiting;                     }
  const std::vector< G4bool >&                       GetWasLimiting()                   const { return fParallelWorldWasLimiting;                     }
  // -- The current and previous step limitation status for the passed parallel world index (obtained with GetParallelWorldIndex(...) above)
  // -- Note that no boundary checks are done on the index passed.
  G4bool                                              GetIsLimiting( G4int worldIndex ) const { return  fParallelWorldIsLimiting[size_t(worldIndex)]; }
  G4bool                                             GetWasLimiting( G4int worldIndex ) const { return fParallelWorldWasLimiting[size_t(worldIndex)]; }
  
  
  // --------------------------------------------------------------
  //                      From process interface
  // --------------------------------------------------------------
  // -- Start/End tracking:
  void StartTracking(G4Track*);
  void   EndTracking();

  // --------------------
  // -- PostStep methods:
  // --------------------
  // -- PostStepGPIL is used to collect up to date volumes in the parallel geometries:
  G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double, G4ForceCondition*);
  // -- PostStepDoIt is not used (never called):
  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step& )
  { return nullptr; }

  // ---------------------------------------------------------------------------
  // -- Along step used for limiting the step on parallel geometries boundaries:
  // ---------------------------------------------------------------------------
  G4double AlongStepGetPhysicalInteractionLength(const G4Track&                track,
						 G4double           previousStepSize,
						 G4double         currentMinimumStep, 
						 G4double&            proposedSafety, 
						 G4GPILSelection*          selection);
  G4VParticleChange*               AlongStepDoIt(const G4Track&                track,
						 const G4Step&                  step);
  
  // ---------------------------
  // -- AtRest methods not used:
  // ---------------------------
  G4double AtRestGetPhysicalInteractionLength(const G4Track&,
					      G4ForceCondition*)
  { return DBL_MAX; }
  
  G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&)
  { return nullptr; }

  
  // --
  virtual void SetProcessManager(const G4ProcessManager*);

    
  
private:
  std::vector< G4VPhysicalVolume* >                      fParallelWorlds;
  std::vector< G4Navigator* >                   fParallelWorldNavigators;
  std::vector< G4int >                    fParallelWorldNavigatorIndeces;
  std::vector< G4double >                         fParallelWorldSafeties;
  std::vector< G4bool >                         fParallelWorldIsLimiting;
  std::vector< G4bool >                        fParallelWorldWasLimiting;
  std::vector< const G4VPhysicalVolume* >                fCurrentVolumes;
  std::vector< const G4VPhysicalVolume* >               fPreviousVolumes;
  G4double                                          fParallelWorldSafety;
  G4bool                                                 fIsTrackingTime;
  G4FieldTrack                                               fFieldTrack;
  G4ParticleChangeForNothing                        fDummyParticleChange;
  G4PathFinder*                                              fPathFinder;
  G4TransportationManager*                        fTransportationManager;
};

#endif
