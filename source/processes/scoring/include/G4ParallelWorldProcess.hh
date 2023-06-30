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
//  G4ParallelWorldProcess.hh
//
//  Description:
//    This procss takes a parallel world and limits a step
//   on the boundaries of volumes in the parallel world.
//    It invokes sensitive detectors assigned in the parallel
//   world.
//    It switches a material (and a region if defined) in the
//   assigned parallel world over the material (and the region)
//   in the mass world.
//
//---------------------------------------------------------------

#ifndef G4ParallelWorldProcess_h
#define G4ParallelWorldProcess_h 1

#include "G4FieldTrack.hh"
#include "G4MultiNavigator.hh"
#include "G4TouchableHandle.hh"
#include "G4VProcess.hh"
#include "globals.hh"

class G4Step;
class G4StepPoint;
class G4Navigator;
class G4TransportationManager;
class G4PathFinder;
class G4VTouchable;
class G4VPhysicalVolume;
class G4ParticleChange;
//------------------------------------------
//
//        G4ParallelWorldProcess class
//
//------------------------------------------

// Class Description:

class G4ParallelWorldProcess : public G4VProcess
{
  public:  // with description
    //------------------------
    // Constructor/Destructor
    //------------------------

    G4ParallelWorldProcess(const G4String& processName = "ParaWorld",
                           G4ProcessType theType = fParallel);
    ~G4ParallelWorldProcess() override;

    //--------------------------------------------------------------
    // Set Paralle World
    //--------------------------------------------------------------

    void SetParallelWorld(G4String parallelWorldName);
    void SetParallelWorld(G4VPhysicalVolume* parallelWorld);

    //--------------------------------------------------------------
    //     Process interface
    //--------------------------------------------------------------

    void StartTracking(G4Track*) override;

    //------------------------------------------------------------------------
    // GetPhysicalInteractionLength() and DoIt() methods for AtRest
    //------------------------------------------------------------------------

    G4double AtRestGetPhysicalInteractionLength(const G4Track&, G4ForceCondition*) override;
    G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&) override;

    //------------------------------------------------------------------------
    // GetPhysicalInteractionLength() and DoIt() methods for AlongStep
    //------------------------------------------------------------------------

    G4double AlongStepGetPhysicalInteractionLength(const G4Track&, G4double, G4double, G4double&,
                                                   G4GPILSelection*) override;
    G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&) override;

    //-----------------------------------------------------------------------
    // GetPhysicalInteractionLength() and DoIt() methods for PostStep
    //-----------------------------------------------------------------------

    G4double PostStepGetPhysicalInteractionLength(const G4Track&, G4double,
                                                  G4ForceCondition*) override;
    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

    //-----------------------------------------------------------------------
    // Flag for material switching
    //-----------------------------------------------------------------------

    inline void SetLayeredMaterialFlag(G4bool flg = true) { layeredMaterialFlag = flg; }
    inline G4bool GetLayeredMaterialFlag() const { return layeredMaterialFlag; }

    //--------------------------------------------------------------------
    // Returns whether a particular particle type requires AtRest process
    //--------------------------------------------------------------------
    G4bool IsAtRestRequired(G4ParticleDefinition*);

    //-----------------------------------------------------------------------
    // Static G4Step object for "Hyper-step"
    //-----------------------------------------------------------------------
    static const G4Step* GetHyperStep();
    static G4int GetHypNavigatorID();

  protected:
    void CopyStep(const G4Step& step);
    void SwitchMaterial(G4StepPoint*);

  protected:
    G4Step* fGhostStep;
    G4StepPoint* fGhostPreStepPoint;
    G4StepPoint* fGhostPostStepPoint;

    G4VParticleChange aDummyParticleChange;
    G4ParticleChange xParticleChange;

    G4TransportationManager* fTransportationManager;
    G4PathFinder* fPathFinder;

    // -------------------------------
    // Navigation in the Ghost World:
    // -------------------------------
    G4String fGhostWorldName;
    G4VPhysicalVolume* fGhostWorld{nullptr};
    G4Navigator* fGhostNavigator{nullptr};
    G4int fNavigatorID{-1};
    G4TouchableHandle fOldGhostTouchable;
    G4TouchableHandle fNewGhostTouchable;
    G4FieldTrack fFieldTrack;
    G4double fGhostSafety{0.};
    G4bool fOnBoundary{false};

    //-----------------------------------------------------------------------
    // Flag for material switching
    //-----------------------------------------------------------------------
    G4bool layeredMaterialFlag{false};

  private:
    static G4ThreadLocal G4Step* fpHyperStep;
    static G4ThreadLocal G4int nParallelWorlds;
    static G4ThreadLocal G4int fNavIDHyp;
    G4int iParallelWorld;
};

#endif
