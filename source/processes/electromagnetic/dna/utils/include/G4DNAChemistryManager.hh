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

#ifndef G4DNACHEMISTRYMANAGER_HH
#define G4DNACHEMISTRYMANAGER_HH

#include "globals.hh"

class G4Track;

enum ElectronicModification
{
    fIonizedMolecule,
    fExcitedMolecule
};

class G4DNAChemistryManager
{
private:
    G4DNAChemistryManager();
    static G4DNAChemistryManager* fInstance;
    bool fActiveChemistry;

public:
    ~G4DNAChemistryManager();
    static G4DNAChemistryManager* Instance();
    static void DeleteInstance();
    inline G4bool IsChemistryActived();
    inline void SetChemistryActivation(G4bool);
    void CreateWaterMolecule(ElectronicModification,
                        G4int /*electronicLevel*/,
                        const G4Track* /*theIncomingTrack*/);
};

inline G4bool G4DNAChemistryManager::IsChemistryActived()
{
    return fActiveChemistry;
}

inline void G4DNAChemistryManager::SetChemistryActivation(G4bool flag)
{
    fActiveChemistry = flag;
}

#endif // G4DNACHEMISTRYMANAGER_HH
