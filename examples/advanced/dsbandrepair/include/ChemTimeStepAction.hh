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
/// \file ChemTimeStepAction.hh
/// \brief Definition of the ChemTimeStepAction class

#ifndef CHEMITACTION_H
#define CHEMITACTION_H

#include "G4UserTimeStepAction.hh"
#include "G4VUserChemistryList.hh"
#include "ChemNtupleManager.hh"

class ChemTimeStepAction : public G4UserTimeStepAction
{
    public:
    ChemTimeStepAction(ChemNtupleManager* histo, TimeStepModel md = fSBS);
    ~ChemTimeStepAction() override = default;
    ChemTimeStepAction(const ChemTimeStepAction& other);
    ChemTimeStepAction& operator=(const ChemTimeStepAction& other);

    void UserReactionAction(const G4Track&a, const G4Track&b,
                                    const std::vector<G4Track*>* products) override;

    void StartProcessing() override;

private:
    ChemNtupleManager* fpHisto{nullptr};

    G4double fReactif1{0.};
    G4double fReactif2{0.};
    G4double fProduct1{0.};
    G4double fProduct2{0.};

    G4double SetFlag(const G4String &val);
};

#endif // ITACTION_H
