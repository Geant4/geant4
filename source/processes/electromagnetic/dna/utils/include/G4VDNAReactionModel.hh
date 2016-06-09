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
// $Id: G4VDNAReactionModel.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
#ifndef G4VReactionModel_
#define G4VReactionModel_

#include "globals.hh"
#include "AddClone_def.hh"

class G4DNAMolecularReactionTable;
class G4Molecule;
class G4Track;

/**
  * G4VDNAReactionModel is an interface used by the G4DNAMolecularReaction process.
  * It defines how the reaction radius should be calculated and whether two molecules
  * can indeed react.
  */

class G4VDNAReactionModel
{
public :
    G4VDNAReactionModel();
    G4VDNAReactionModel(const G4VDNAReactionModel&);
    virtual ~G4VDNAReactionModel();

    /** This macro is defined in AddClone_def **/
    G4IT_TO_BE_CLONED(G4VDNAReactionModel)

    virtual void Initialise(const G4Molecule*, const G4Track&) {;}
    virtual void InitialiseToPrint(const G4Molecule*) = 0 ;
    virtual G4double GetReactionRadius(const G4Molecule*, const G4Molecule*) = 0;
    virtual G4double GetReactionRadius(const int) = 0;
    virtual G4bool FindReaction(const G4Track&, const G4Track&,
                                const G4double /*reactionRadius*/,
                                G4double& /*separationDistance*/,  // To be calculated
                                const G4bool /*hasReachedUserTimeLimit*/) = 0;

    inline void SetReactionTable(const G4DNAMolecularReactionTable*);
    inline const G4DNAMolecularReactionTable* GetReactionTable();

protected :
    G4VDNAReactionModel& operator=(const G4VDNAReactionModel&);
    const G4DNAMolecularReactionTable* fReactionTable ;
};

inline void G4VDNAReactionModel::SetReactionTable(const G4DNAMolecularReactionTable* table)
{
    fReactionTable = table ;
}

inline const G4DNAMolecularReactionTable* G4VDNAReactionModel::GetReactionTable()
{
    return fReactionTable ;
}
#endif
