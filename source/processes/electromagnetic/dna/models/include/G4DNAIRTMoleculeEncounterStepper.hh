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
/*
 * G4DNAIRTMoleculeEncounterStepper.hh
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

#pragma once

#include "G4VITTimeStepComputer.hh"
#include "G4KDTreeResult.hh"
#include "G4ITReaction.hh"
#include "G4ITTrackHolder.hh"

class G4VDNAReactionModel;
class G4DNAMolecularReactionTable;
class G4MolecularConfiguration;

class G4Molecule;

/**
 * Given a molecule G4DNAIRTMoleculeEncounterStepper will calculate for its possible reactants
 * what will be the minimum encounter time and the associated molecules.*
 *
 * This model includes dynamical time steps as explained in
 * "Computer-Aided Stochastic Modeling of the Radiolysis of Liquid Water",
 * V. Michalik, M. Begusová, E. A. Bigildeev,
 * Radiation Research, Vol. 149, No. 3 (Mar., 1998), pp. 224-236
 *
 */

class G4DNAIRTMoleculeEncounterStepper : public G4VITTimeStepComputer
{
public:
    G4DNAIRTMoleculeEncounterStepper();
    virtual ~G4DNAIRTMoleculeEncounterStepper();
    G4DNAIRTMoleculeEncounterStepper(const G4DNAIRTMoleculeEncounterStepper&) = delete;
    G4DNAIRTMoleculeEncounterStepper& operator=(const G4DNAIRTMoleculeEncounterStepper&) = delete;

    virtual void Prepare();
    virtual G4double CalculateStep(const G4Track&, const G4double&);
    virtual G4double CalculateMinTimeStep(G4double, G4double);

    void SetReactionModel(G4VDNAReactionModel*);
    G4VDNAReactionModel* GetReactionModel();

    void SetVerbose(int);
    // Final time returned when reaction is available in the reaction table = 1
    // All details = 2

private:
    void InitializeForNewTrack();

    class Utils;
    void CheckAndRecordResults(const Utils&,
#ifdef G4VERBOSE
                               const G4double reactionRange,
#endif
                               G4KDTreeResultHandle&);

    G4bool fHasAlreadyReachedNullTime;

    const G4DNAMolecularReactionTable*& fMolecularReactionTable;
    G4VDNAReactionModel* fReactionModel;
    G4ITReactionSet* fReactionSet;
    G4ITTrackHolder* fpTrackContainer;
    G4int fVerbose;

    class Utils
    {
    public:
        Utils(const G4Track& tA, const G4MolecularConfiguration* mB);
        ~Utils() = default;

        G4double GetConstant() const
        {
            return fConstant;
        }

        const G4Track& fpTrackA;
        const G4MolecularConfiguration* fpMoleculeB;
        const G4Molecule* fpMoleculeA;
        G4double fDA;
        G4double fDB;
        G4double fConstant;
    };
};
