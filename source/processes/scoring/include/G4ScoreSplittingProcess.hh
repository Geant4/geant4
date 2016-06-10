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
// $Id: G4ScoreSplittingProcess.hh 68733 2013-04-05 09:45:28Z gcosmo $
//
// 
//---------------------------------------------------------------
//
//  G4ScoreSplittingProcess.hh
//
//  Description:
//    This process is used to split the length and energy 
//   of a step in a regular structure into sub-steps, and to 
//   call the scorers for each sub-volume.
//    It invokes sensitive detectors assigned in the *mass*
//   world.
//
//  Design and first implementation: J. Apostolakis / M.Asai 2010
//---------------------------------------------------------------


#ifndef G4ScoreSplittingProcess_h
#define G4ScoreSplittingProcess_h 1

#include "globals.hh"
class G4Step;
class G4Navigator;
class G4TransportationManager;
class G4PathFinder;
class G4VTouchable;
class G4VPhysicalVolume;
class G4ParticleChange;
class G4EnergySplitter; 

#include "G4VProcess.hh"
#include "G4FieldTrack.hh"
#include "G4TouchableHandle.hh"
class G4TouchableHistory;
// #include "G4TouchableHistory.hh"

//------------------------------------------
//
//        G4ScoreSplittingProcess class
//
//------------------------------------------


// Class Description:

class G4ScoreSplittingProcess : public G4VProcess
{
public: // with description

  //------------------------
  // Constructor/Destructor
  //------------------------
  
  G4ScoreSplittingProcess(const G4String& processName = "ScoreSplittingProc",
				 G4ProcessType theType = fParameterisation);
  virtual ~G4ScoreSplittingProcess();
  
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
  G4TouchableHistory* CreateTouchableForSubStep( G4int newVoxelNum, G4ThreeVector newPosition ); 

private:
  void CopyStepStart(const G4Step & step);

  G4Step * fSplitStep;
  G4StepPoint *fSplitPreStepPoint;
  G4StepPoint *fSplitPostStepPoint;

  G4VParticleChange dummyParticleChange;
  G4ParticleChange xParticleChange;

  // G4TransportationManager* fTransportationManager;
  // G4PathFinder*        fPathFinder;

  // -------------------------------
  // Touchables for the Split Step
  // -------------------------------
  G4TouchableHandle    fOldTouchableH;
  G4TouchableHandle    fNewTouchableH;

  // Memory of Touchables of full step
  G4TouchableHandle    fInitialTouchableH;
  G4TouchableHandle    fFinalTouchableH;

  G4EnergySplitter     *fpEnergySplitter; 

  // ******************************************************
  //  For TESTS:
  // ******************************************************
public:
  void Verbose(const G4Step&) const;
};

#endif
