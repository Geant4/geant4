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
// $Id: $
//---------------------------------------------------------------
//
// G4ChannelingOptrChangeCrossSection
// Based On GB01BOptrChangeCrossSection
// Class Description:
//        A G4VBiasingOperator concrete implementation example to
//    illustrate how to bias physics processes cross-section for
//    one particle type.
//        The G4VBiasingOperation G4BOptnChangeCrossSection is
//    selected by this operator, and is sent to each process
//    calling the operator.
//        A simple constant bias to the cross-section is applied,
//    but more sophisticated changes can be applied.
//
//---------------------------------------------------------------
//

#ifndef G4ChannelingOptrChangeCrossSection_hh
#define G4ChannelingOptrChangeCrossSection_hh 1

#include "G4VBiasingOperator.hh"
class G4BOptnChangeCrossSection;
class G4ParticleDefinition;
#include <map>
#include <unordered_map>

enum G4ChannelingDensityRatio{
    fDensityRatioNotDefined = -1,
    fDensityRatioNone = 0,
    fDensityRatioNuDElD = 1,
    fDensityRatioNuD = 2,
    fDensityRatioElD = 3
};

class G4ChannelingOptrChangeCrossSection : public G4VBiasingOperator {
public:
    // ------------------------------------------------------------
    // -- Constructor: takes the name of the particle type to bias:
    // ------------------------------------------------------------
    G4ChannelingOptrChangeCrossSection(G4String particleToBias, G4String name = "ChannelingChangeXS");
    virtual ~G4ChannelingOptrChangeCrossSection();
    virtual void StartRun();
    
private:
    G4int fChannelingID;
    
private:
    virtual G4VBiasingOperation*
    ProposeOccurenceBiasingOperation(const G4Track*                            track,
                                     const G4BiasingProcessInterface* callingProcess);
    virtual G4VBiasingOperation*
    ProposeFinalStateBiasingOperation(const G4Track*, const G4BiasingProcessInterface*)
    {return 0;}
    virtual G4VBiasingOperation*
    ProposeNonPhysicsBiasingOperation(const G4Track*, const G4BiasingProcessInterface*)
    {return 0;}
    
private:
    using G4VBiasingOperator::OperationApplied;
    virtual void OperationApplied( const G4BiasingProcessInterface*                callingProcess,
                                  G4BiasingAppliedCase                               biasingCase,
                                  G4VBiasingOperation*                 occurenceOperationApplied,
                                  G4double                         weightForOccurenceInteraction,
                                  G4VBiasingOperation*                finalStateOperationApplied,
                                  const G4VParticleChange*                particleChangeProduced );
    
private:
    std::map< const G4BiasingProcessInterface*,
    G4BOptnChangeCrossSection*       > fChangeCrossSectionOperations;
    G4bool                                  fSetup;
    const G4ParticleDefinition*    fParticleToBias;
    
private:
    std::unordered_map<std::string,G4ChannelingDensityRatio> fProcessToDensity;
};

#endif
