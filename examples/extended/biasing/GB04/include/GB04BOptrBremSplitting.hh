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
///
/// \file GB04BOptrBremSplitting.hh
/// \brief Definition of the GB04BOptrBremSplitting class
#ifndef GB04BOptrBremSplitting_hh
#define GB04BOptrBremSplitting_hh 1

#include "G4VBiasingOperator.hh"
class GB04BOptnBremSplitting;
class G4GenericMessenger;

class GB04BOptrBremSplitting : public G4VBiasingOperator {
public:
  GB04BOptrBremSplitting();
  virtual ~GB04BOptrBremSplitting() {}
  
public:
  // -------------------------
  // Optional from base class:
  // -------------------------
  // -- Call at run start:
  virtual void      StartRun();
  // -- Call at each track starting:
  virtual void StartTracking( const G4Track* track );

private:
  // -----------------------------
  // -- Mandatory from base class:
  // -----------------------------
  // -- Unused:
  virtual G4VBiasingOperation*
  ProposeNonPhysicsBiasingOperation(const G4Track* /* track */,
                                    const G4BiasingProcessInterface* /* callingProcess */)
  { return 0; }
  virtual G4VBiasingOperation* 
  ProposeOccurenceBiasingOperation (const G4Track* /* track */,
                                    const G4BiasingProcessInterface* /* callingProcess */)
  { return 0; }
  // -- Used:
  virtual G4VBiasingOperation*
  ProposeFinalStateBiasingOperation(const G4Track* track,
                                    const G4BiasingProcessInterface* callingProcess);
  
private:
  // -- Avoid compiler complaining for (wrong) method shadowing,
  // -- this is because other virtual method with same name exists.
  using G4VBiasingOperator::OperationApplied;

private:
  GB04BOptnBremSplitting* fBremSplittingOperation;
  G4int                          fSplittingFactor;
  G4bool                         fBiasPrimaryOnly;
  G4bool                            fBiasOnlyOnce;
  G4int                            fNInteractions;
  // Messengers to change the 
  G4GenericMessenger*  fSplittingFactorMessenger;
  G4GenericMessenger*  fBiasPrimaryOnlyMessenger;
  G4GenericMessenger*     fBiasOnlyOnceMessenger;
  
  
  
};

#endif
