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
/// \file GB07/include/GB07BOptrLeadingParticle.hh
/// \brief Definition of the GB07BOptrLeadingParticle class
//
#ifndef GB07BOptrLeadingParticle_hh
#define GB07BOptrLeadingParticle_hh 1

#include "G4VBiasingOperator.hh"
class G4BOptnLeadingParticle;
class G4VProcess;

#include <map>

class GB07BOptrLeadingParticle : public G4VBiasingOperator {
public:
  GB07BOptrLeadingParticle( G4String operatorName = "LeadingParticleBiasingOperator");
  virtual ~GB07BOptrLeadingParticle();
  
private:
  // -----------------------------
  // -- Mandatory from base class:
  // -----------------------------
  // -- Unsused:
  virtual G4VBiasingOperation*
  ProposeNonPhysicsBiasingOperation(const G4Track*,
                                    const G4BiasingProcessInterface*) final
  { return nullptr; }
  // -- Unused:
  virtual G4VBiasingOperation* 
  ProposeOccurenceBiasingOperation (const G4Track*,
                                    const G4BiasingProcessInterface*) final
  { return nullptr; }
  // -- Used:
  // -- Will return the biasing operation at the final state generation stage
  virtual G4VBiasingOperation*
  ProposeFinalStateBiasingOperation(const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess) final;

public:
  // ------------------------------------
  // -- Optional methods from base class:
  // ------------------------------------
  virtual void      StartRun()                       final;
  virtual void StartTracking( const G4Track* track ) final;

private:
  // -- The leading particle biasing operation that will actually
  // -- trim the final state generation:
  G4BOptnLeadingParticle* fLeadingParticleBiasingOperation;
  // -- Two-particle final states, they will be biased with a
  // -- "further killing probability" for particles accompanying
  // -- the leading particle, we remember what processes are concerned:
  const G4VProcess*       fAnnihilation;
  const G4VProcess*         fConversion;
  const G4VProcess*              fDecay;
  const G4VProcess* fTwoParticleProcess;
  
};

#endif
