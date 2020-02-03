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
/// \file BiasingOperation.hh
/// \brief Definition of the BiasingOperation class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef BiasingOperation_hh
#define BiasingOperation_hh 1

#include "G4VBiasingOperation.hh"

class G4ProtonInelasticProcess;
class G4NeutronInelasticProcess;
class G4PionPlusInelasticProcess;
class G4PionMinusInelasticProcess;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BiasingOperation : public G4VBiasingOperation {
  // The biasing operation implemented in this class is indeed a "trick" to 
  // use FTFP+INCLXX instead of FTFP+BERT for determining the final-state of
  // proton, neutron, pion+, pion- inelastic interactions happening in one
  // particular logical volume, Tracking_region, where the biasing is applied.
  public:
    BiasingOperation( G4String name );
    virtual ~BiasingOperation();
    virtual G4VParticleChange* ApplyFinalStateBiasing( const G4BiasingProcessInterface*, 
                                                       const G4Track*, const G4Step*, G4bool& );
    // Unused :
    virtual const G4VBiasingInteractionLaw* 
      ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface*, 
                                             G4ForceCondition& ) { return 0; }
    virtual G4double 
      DistanceToApplyOperation( const G4Track*, G4double, G4ForceCondition* ) { return DBL_MAX; }
    virtual G4VParticleChange*
      GenerateBiasingFinalState( const G4Track*, const G4Step* ) { return 0; }

  private:
    G4ProtonInelasticProcess*    fProtonInelasticProcess;
    G4NeutronInelasticProcess*   fNeutronInelasticProcess;
    G4PionPlusInelasticProcess*  fPionPlusInelasticProcess;
    G4PionMinusInelasticProcess* fPionMinusInelasticProcess;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

