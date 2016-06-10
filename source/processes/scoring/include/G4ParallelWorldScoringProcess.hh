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
// $Id: G4ParallelWorldScoringProcess.hh 68733 2013-04-05 09:45:28Z gcosmo $
//
// 
//---------------------------------------------------------------
//
//  G4ParallelWorldScoringProcess.hh
//
//  Description:
//    This procss takes a parallel world and limits a step
//   on the boundaries of volumes in the parallel world.
//    It invokes sensitive detectors assigned in the parallel
//   world.
//
//---------------------------------------------------------------


#ifndef G4ParallelWorldScoringProcess_h
#define G4ParallelWorldScoringProcess_h 1

#include "globals.hh"
class G4Step;
class G4Navigator;
class G4TransportationManager;
class G4PathFinder;
class G4VTouchable;
class G4VPhysicalVolume;
class G4ParticleChange;
class G4ParticleDefinition;
#include "G4VProcess.hh"
#include "G4FieldTrack.hh"
#include "G4TouchableHandle.hh"

//------------------------------------------
//
//        G4ParallelWorldScoringProcess class
//
//------------------------------------------


// Class Description:

class G4ParallelWorldScoringProcess : public G4VProcess
{
public: // with description

  //------------------------
  // Constructor/Destructor
  //------------------------
  
  G4ParallelWorldScoringProcess(const G4String& processName = "ParaWorldScore",
				 G4ProcessType theType = fParameterisation);
  virtual ~G4ParallelWorldScoringProcess();
  
  //--------------------------------------------------------------
  // Set Paralle World
  //--------------------------------------------------------------

  void SetParallelWorld(G4String parallelWorldName);
  void SetParallelWorld(G4VPhysicalVolume* parallelWorld);
  G4bool IsAtRestRequired(G4ParticleDefinition* partDef);

  //--------------------------------------------------------------
  //     Process interface
  //--------------------------------------------------------------

  void StartTracking(G4Track*);
  
  //------------------------------------------------------------------------
  // GetPhysicalInteractionLength() and DoIt() methods for AtRest 
  //------------------------------------------------------------------------
  
  G4double AtRestGetPhysicalInteractionLength(
					      const G4Track& ,
					      G4ForceCondition* 
					      );

  G4VParticleChange* AtRestDoIt(
			       const G4Track& ,
			       const G4Step&
			       );

  //------------------------------------------------------------------------
  // GetPhysicalInteractionLength() and DoIt() methods for AlongStep 
  //------------------------------------------------------------------------
  
  G4double AlongStepGetPhysicalInteractionLength(
						 const G4Track&,
						 G4double  ,
						 G4double  ,
						 G4double&,
						 G4GPILSelection*
						 );

  G4VParticleChange* AlongStepDoIt(
				  const G4Track& ,
				  const G4Step& 
				  );

  //-----------------------------------------------------------------------
  // GetPhysicalInteractionLength() and DoIt() methods for PostStep
  //-----------------------------------------------------------------------
  
  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
						G4double   previousStepSize,
						G4ForceCondition* condition);
  
  G4VParticleChange* PostStepDoIt(const G4Track& ,const G4Step& );

private:
  void CopyStep(const G4Step & step);

  G4Step * fGhostStep;
  G4StepPoint * fGhostPreStepPoint;
  G4StepPoint * fGhostPostStepPoint;

  G4VParticleChange aDummyParticleChange;
  G4ParticleChange xParticleChange;

  G4TransportationManager* fTransportationManager;
  G4PathFinder*        fPathFinder;

  // -------------------------------
  // Navigation in the Ghost World:
  // -------------------------------
  G4String             fGhostWorldName;
  G4VPhysicalVolume*   fGhostWorld;
  G4Navigator*         fGhostNavigator;
  G4int                fNavigatorID;
  G4TouchableHandle    fOldGhostTouchable;
  G4TouchableHandle    fNewGhostTouchable;
  G4FieldTrack         fFieldTrack;
  G4double             fGhostSafety;
  G4bool               fOnBoundary;

  // ******************************************************
  // ******************************************************
  //
  //  For TESTS:
  //
  // ******************************************************
  // ******************************************************
public:
  void Verbose(const G4Step&) const;
};

#endif
