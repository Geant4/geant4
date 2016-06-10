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
                                    const G4BiasingProcessInterface* callingProcess);
  virtual G4VBiasingOperation* 
  ProposeOccurenceBiasingOperation (const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess);
  virtual G4VBiasingOperation*
  ProposeFinalStateBiasingOperation(const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess);
  
private:
  // -- Avoid compiler complaining for (wrong) method shadowing,
  // -- this is because other virtual method with same name exists.
  using G4VBiasingOperator::OperationApplied;
  
  // -- Optionnal, from base class:
  // -- Here, nedded implemention to forward calls from process to the underneath
  // -- G4BOptrForceCollision biasing operator
  virtual void OperationApplied( const G4BiasingProcessInterface*         callingProcess,
                                 G4BiasingAppliedCase                        biasingCase,
                                 G4VBiasingOperation*                   operationApplied,
                                 const G4VParticleChange*         particleChangeProduced );
  void ExitBiasing( const G4Track*, const G4BiasingProcessInterface* );
  
public:
  virtual void StartTracking( const G4Track* track );
  
private:
  std::map < const G4ParticleDefinition*, G4BOptrForceCollision* > fBOptrForParticle;
  std::vector < const G4ParticleDefinition* > fParticlesToBias;
  G4BOptrForceCollision* fCurrentOperator;
};

#endif
