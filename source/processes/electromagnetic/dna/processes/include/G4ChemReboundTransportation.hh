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

#ifndef TRANSPORTATION_H
#define TRANSPORTATION_H

#include "G4DNABoundingBox.hh"
#include "G4ITReaction.hh"
#include "G4ITTransportation.hh"
#include "G4MolecularConfiguration.hh"

class G4Molecule;
class G4VUserBrownianAction;

class G4ChemReboundTransportation : public G4ITTransportation
{
    using MolConf = const G4MolecularConfiguration*;

public:
    explicit G4ChemReboundTransportation(const G4String& aName = "ChemReboundTransportation",
                                       const G4DNABoundingBox* = nullptr, G4int verbosityLevel = 0);
    ~G4ChemReboundTransportation() override = default;
    G4ChemReboundTransportation(const G4ChemReboundTransportation&) = delete;
    G4ChemReboundTransportation& operator=(const G4ChemReboundTransportation&) = delete;
    void BuildPhysicsTable(const G4ParticleDefinition&) override;
    void StartTracking(G4Track* aTrack) override;
    void ComputeStep(const G4Track&, const G4Step&, G4double, G4double&) override;
    G4double AlongStepGetPhysicalInteractionLength(const G4Track& /*track*/,
                                                   G4double /*previousStepSize*/,
                                                   G4double /*currentMinimumStep*/,
                                                   G4double& /*currentSafety*/,
                                                   G4GPILSelection* /*selection*/) override;

    G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step&) override;
    G4VParticleChange* AlongStepDoIt(const G4Track& track, const G4Step&) override;

    inline void SetBoundary(const G4DNABoundingBox*);
    static G4double calculateNextCoordinate(G4double nextPos, G4double high, G4double low);
    G4double GetTimeToBoundary(const G4Track& track);

protected:
    //________________________________________________________________
    // Process information
    struct G4ITBrownianState : public G4ITTransportationState
    {
    public:
        G4ITBrownianState();
        ~G4ITBrownianState() override = default;
        G4String GetType() override { return "Rebound"; }

        G4bool fTimeStepReachedLimit;
        G4double fRandomNumber;
    };
    G4ITReactionSet* fReactionSet = G4ITReactionSet::Instance();

    G4Material* fNistWater = nullptr;
    G4double fMaximumTimeStep = 0;
    G4double fInternalMinTimeStep;
    const G4DNABoundingBox* fpBoundingBox = nullptr;
    G4ThreeVector BouncingAction(const G4ThreeVector& nextPosition);
    G4double calculateDistanceFromTimeStep(MolConf mol, G4double timeStep);
};

inline void G4ChemReboundTransportation::SetBoundary(const G4DNABoundingBox* pBounding)
{
  fpBoundingBox = pBounding;
}

#endif
