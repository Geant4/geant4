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
// $Id: G4DNAMoleculeEncounterStepper.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4MOLECULEENCOUNTERSTEPPER_H
#define G4MOLECULEENCOUNTERSTEPPER_H

#include "G4VITTimeStepper.hh"
#include "G4ITManager.hh"

class G4VDNAReactionModel;
class G4DNAMolecularReactionTable;

class G4Molecule;

/**
  * Given a molecule G4DNAMoleculeEncounterStepper will calculate for its possible reactants
  * what will be the minimum encounter time and the associated molecules.*
  *
  * This model includes dynamical time steps as explained in
  * "Computer-Aided Stochastic Modeling of the Radiolysis of Liquid Water",
  * V. Michalik, M. Begusová, E. A. Bigildeev,
  * Radiation Research, Vol. 149, No. 3 (Mar., 1998), pp. 224-236
  *
  */

class G4DNAMoleculeEncounterStepper : public G4VITTimeStepper
{
public:
    G4DNAMoleculeEncounterStepper();
    virtual ~G4DNAMoleculeEncounterStepper();
    G4DNAMoleculeEncounterStepper(const G4DNAMoleculeEncounterStepper&);
    G4IT_ADD_CLONE(G4VITTimeStepper,G4DNAMoleculeEncounterStepper)

    virtual void Prepare();
//    virtual void PrepareForAllProcessors();
    virtual G4double CalculateStep(const G4Track&, const G4double&);

    inline void SetReactionModel(G4VDNAReactionModel*);
    inline G4VDNAReactionModel* GetReactionModel();

    inline void SetVerbose(int);
    // Final time returned when reaction is avalaible in the reaction table = 1
    // All details = 2

private:

    void RetrieveResults(const G4Track&, const G4Molecule*, const G4Molecule*, const G4double /*reactionRange*/,
                         G4KDTreeResultHandle&, G4bool iterate = true);

    G4bool fHasAlreadyReachedNullTime;

    G4DNAMoleculeEncounterStepper& operator=(const G4DNAMoleculeEncounterStepper&);
    const G4DNAMolecularReactionTable*& fMolecularReactionTable ;
    G4VDNAReactionModel* fReactionModel;
    G4int fVerbose ;
};

inline void G4DNAMoleculeEncounterStepper::SetReactionModel(G4VDNAReactionModel* reactionModel)
{
    fReactionModel = reactionModel ;
}

inline G4VDNAReactionModel* G4DNAMoleculeEncounterStepper::GetReactionModel()
{
    return fReactionModel;
}

inline void G4DNAMoleculeEncounterStepper::SetVerbose(int flag)
{
    fVerbose = flag;
}

#endif // G4MOLECULEENCOUNTERSTEPPER_H
