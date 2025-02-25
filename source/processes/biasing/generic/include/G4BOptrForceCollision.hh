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
// G4BOptrForceCollision
//
// Class Description:
//
// A G4VBiasingOperator that implements a "force collision" a la
// MCNP. This is meant for neutral particles.
// When the track enters the volume, it is cloned. One copy makes
// a forced free flight up to the volume exit. The other copy makes
// a forced collision inside the volume.
//
// Author: Marc Verderi, November 2013.
// --------------------------------------------------------------------
#ifndef G4BOptrForceCollision_hh
#define G4BOptrForceCollision_hh 1

#include "G4VBiasingOperator.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include <map>

class G4BOptnForceFreeFlight;
class G4BOptnForceCommonTruncatedExp;
class G4BOptnCloning;
class G4VProcess;
class G4BiasingProcessInterface;
class G4ParticleDefinition;
class G4BOptrForceCollisionTrackData;

class G4BOptrForceCollision : public G4VBiasingOperator
{
  public:

    G4BOptrForceCollision(const G4String& particleToForce,
                          const G4String& name="ForceCollision");
    G4BOptrForceCollision(const G4ParticleDefinition* particleToForce,
                          const G4String& name="ForceCollision");
    ~G4BOptrForceCollision();
  
    virtual void Configure() final;
    virtual void ConfigureForWorker() final;
    virtual void StartRun() final;
    virtual void StartTracking( const G4Track* track ) final;
    virtual void ExitBiasing( const G4Track*, const G4BiasingProcessInterface* ) final {};
    virtual void EndTracking() final;

    // -- operation applied:
    void OperationApplied( const G4BiasingProcessInterface* callingProcess,
                                 G4BiasingAppliedCase biasingCase,
                                 G4VBiasingOperation* operationApplied,
                           const G4VParticleChange* particleChangeProduced ) final;
    void OperationApplied( const G4BiasingProcessInterface* callingProcess,
                                 G4BiasingAppliedCase biasingCase,
                                 G4VBiasingOperation* occurenceOperationApplied,
                                 G4double weightForOccurenceInteraction,
                                 G4VBiasingOperation* finalStateOperationApplied,
                           const G4VParticleChange* particleChangeProduced ) final;

  private:

    // -- Mandatory from base class :
    virtual G4VBiasingOperation*
    ProposeNonPhysicsBiasingOperation(const G4Track* track,
                         const G4BiasingProcessInterface* callingProcess) final;
    virtual G4VBiasingOperation*
    ProposeOccurenceBiasingOperation(const G4Track* track,
                         const G4BiasingProcessInterface* callingProcess) final;
    virtual G4VBiasingOperation*
    ProposeFinalStateBiasingOperation(const G4Track* track,
                         const G4BiasingProcessInterface* callingProcess) final;

  private:

    G4int fForceCollisionModelID = 0;
    const G4Track* fCurrentTrack = nullptr;
    G4BOptrForceCollisionTrackData* fCurrentTrackData = nullptr;
    std::map< const G4BiasingProcessInterface*, G4BOptnForceFreeFlight* > fFreeFlightOperations;
    G4BOptnForceCommonTruncatedExp* fSharedForceInteractionOperation = nullptr;
    G4BOptnCloning* fCloningOperation = nullptr;
    G4double fInitialTrackWeight = -1.0;
    G4bool fSetup = true;
    const G4ParticleDefinition* fParticleToBias = nullptr;
};

#endif
