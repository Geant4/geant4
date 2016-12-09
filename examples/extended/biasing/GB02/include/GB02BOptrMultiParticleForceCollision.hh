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
/// \file GB02/include/GB02BOptrMultiParticleForceCollision.hh
/// \brief Definition of the GB02BOptrMultiParticleForceCollision class
//
#ifndef GB02BOptrMultiParticleForceCollision_hh
#define GB02BOptrMultiParticleForceCollision_hh 1

#include "G4VBiasingOperator.hh"
class G4BOptrForceCollision;
class G4ParticleDefinition;

#include <map>

class GB02BOptrMultiParticleForceCollision : public G4VBiasingOperator {
public:
  GB02BOptrMultiParticleForceCollision();
  virtual ~GB02BOptrMultiParticleForceCollision() {}

  // --------------------------
  // -- Specific to this class:
  // --------------------------
  // -- Declare particles to be biased:
  void AddParticle( G4String particleName );
  
private:
  // -----------------------------
  // -- Mandatory from base class:
  // -----------------------------
  virtual G4VBiasingOperation*
  ProposeNonPhysicsBiasingOperation(const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess) final;
  virtual G4VBiasingOperation* 
  ProposeOccurenceBiasingOperation (const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess) final;
  virtual G4VBiasingOperation*
  ProposeFinalStateBiasingOperation(const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess) final;
  
private:
  // -------------------------------
  // -- Optionnal, from base class :
  // -------------------------------
  // -- Here, these are needed implementions to forward calls to the underneath
  // -- G4BOptrForceCollision biasing operator
  void OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                         G4BiasingAppliedCase                        biasingCase,
                         G4VBiasingOperation*                   operationApplied,
                         const G4VParticleChange*         particleChangeProduced ) final;
  void OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                         G4BiasingAppliedCase                        biasingCase,
                         G4VBiasingOperation*          occurenceOperationApplied,
                         G4double                  weightForOccurenceInteraction,
                         G4VBiasingOperation*         finalStateOperationApplied, 
                         const G4VParticleChange*         particleChangeProduced ) final;
  
  void ExitBiasing( const G4Track*, const G4BiasingProcessInterface* ) final;
  
public:
  virtual void StartTracking( const G4Track* track ) final;
  
private:
  std::map < const G4ParticleDefinition*, G4BOptrForceCollision* > fBOptrForParticle;
  std::vector < const G4ParticleDefinition* > fParticlesToBias;
  G4BOptrForceCollision* fCurrentOperator;
};

#endif
