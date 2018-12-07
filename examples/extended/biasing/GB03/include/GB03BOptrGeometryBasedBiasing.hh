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
/// \file GB03BOptrGeometryBasedBiasing.hh
/// \brief Definition of the GB03BOptrGeometryBasedBiasing class

#ifndef GB03BOptrGeometryBasedBiasing_h
#define GB03BOptrGeometryBasedBiasing_h 1

#include "G4VBiasingOperator.hh"
#include "GB03BOptnSplitOrKillOnBoundary.hh"
class G4GenericMessenger;

/// Biasing operator class.

class GB03BOptrGeometryBasedBiasing : public G4VBiasingOperator
{
  public:
    GB03BOptrGeometryBasedBiasing();
    virtual ~GB03BOptrGeometryBasedBiasing();

public:
  // ------------------------------
  // Method added for this example:
  // ------------------------------
  GB03BOptnSplitOrKillOnBoundary* GetSplitAndKillOperation() const
  { return fSplitAndKillOperation; }

  // -------------------------
  // Optional from base class:
  // -------------------------
  void StartRun();
  
private:
  // --------------------------
  // Mandatory from base class:
  // --------------------------
  // Used for splitting/killing:
  virtual G4VBiasingOperation*
  ProposeNonPhysicsBiasingOperation( const G4Track* track,
                                     const G4BiasingProcessInterface* callingProcess );

  // Not used here:
  virtual G4VBiasingOperation* 
  ProposeOccurenceBiasingOperation( const G4Track*,
                                    const G4BiasingProcessInterface* ) { return 0; }
  // Not used here:
  virtual G4VBiasingOperation*
  ProposeFinalStateBiasingOperation( const G4Track*,
                                     const G4BiasingProcessInterface* ) { return 0; }

private:
  GB03BOptnSplitOrKillOnBoundary* fSplitAndKillOperation;
  G4int    fSplittingFactor;
  G4double fApplyProbability;
  // Messengers to change the 
  G4GenericMessenger*  fSplittingFactorMessenger;
  G4GenericMessenger* fApplyProbabilityMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
