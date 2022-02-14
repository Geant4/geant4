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
/*
 * G4DNAElectronHoleRecombination.hh
 *
 *  Created on: Jun 17, 2015
 *      Author: mkaramit
 */

#pragma once
#include <G4VITRestDiscreteProcess.hh>

class G4DNAElectronHoleRecombination : public G4VITRestDiscreteProcess
{
public:
    G4DNAElectronHoleRecombination();
    virtual ~G4DNAElectronHoleRecombination();
    void Create();

    void StartTracking(G4Track*) override;

    G4bool IsApplicable(const G4ParticleDefinition&) override;

    void BuildPhysicsTable(const G4ParticleDefinition&) override
    {}

    ////////////////////////////
    // DoIt    /////////////////
    ///////////////////////////
    G4VParticleChange* AtRestDoIt(const G4Track& /*track*/,
                                          const G4Step& /*stepData*/) override;
    //  A virtual base class function that has to be overridden
    //  by any subclass. The DoIt method actually performs the
    //  physics process and determines either momentum change
    //  of the production of secondaries etc.
    //    arguments
    //      const G4Track&    track:
    //        reference to the current G4Track information
    //      const G4Step&     stepData:
    //        reference to the current G4Step information

    G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;


    //////////////////////////
    // GPIL    //////////////
    /////////////////////////

    G4double GetMeanFreePath(const G4Track& aTrack,
                                    G4double previousStepSize,
                                    G4ForceCondition* condition)override;
    //  Calculates from the macroscopic cross section a mean
    //  free path, the value is returned in units of distance.

    G4double GetMeanLifeTime(const G4Track& aTrack,
                                    G4ForceCondition* condition) override;
    //  Calculates the mean life-time (i.e. for decays) of the
    //  particle at rest due to the occurrence of the given process,
    //  or converts the probability of interaction (i.e. for
    //  annihilation) into the life-time of the particle for the
    //  occurence of the given process.

private:
    struct ReactantInfo
    {
        G4Track* fElectron;
        G4double fDistance;
        G4double fProbability;
    };

    struct State : public G4ProcessStateBase<G4DNAElectronHoleRecombination>
    {
        std::vector<ReactantInfo> fReactants; // distanceSqr;
        G4double fSampleProba;
    };

    G4bool FindReactant(const G4Track& track);
    void MakeReaction(const G4Track& track);

    const std::vector<double>* fpMoleculeDensity;
    G4ParticleChange fParticleChange;
    G4bool fIsInitialized;
    std::map<int, std::pair<double, G4double> > fOnsagerRadiusPerMaterial;
};
