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
/// \file BiasingOperator.hh
/// \brief Definition of the BiasingOperator class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BiasingOperator_hh
#define BiasingOperator_hh 1

#include "G4VBiasingOperator.hh"
#include <vector>

class G4ParticleDefinition;
class BiasingOperation;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BiasingOperator : public G4VBiasingOperator {
  // When a proton, or neutron, or pion+ or pion- inelastic process occurs
  // (naturally, without any biasing) in the logical volume(s) where this
  // biasing operator has been attached to, this class uses the biasing "trick"
  // of calling FTFP+INCLXX instead of FTFP+BERT for determining the final-state.
  // Note that the weights of the produced secondaries are left to their default values, 1.0.
  public:
    BiasingOperator();
    virtual ~BiasingOperator() {}
    void AddParticle( G4String particleName );
    virtual G4VBiasingOperation* ProposeFinalStateBiasingOperation( const G4Track* track,
                                        const G4BiasingProcessInterface* callingProcess ) final;
    // Not used:  
    virtual G4VBiasingOperation* ProposeNonPhysicsBiasingOperation( const G4Track*,
                                    const G4BiasingProcessInterface* ) { return 0; }
    virtual G4VBiasingOperation* ProposeOccurenceBiasingOperation ( const G4Track*,
                                    const G4BiasingProcessInterface* ) { return 0; }
  private:
    std::vector< const G4ParticleDefinition* > fParticlesToBias;
    BiasingOperation* fBiasingOperation;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

