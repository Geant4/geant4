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


#ifndef G4DNAMakeReaction_hh
#define G4DNAMakeReaction_hh 1

#include "G4VITReactionProcess.hh"
class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;
class G4ITReactionSet;
class G4VITTimeStepComputer;

class G4DNAMakeReaction : public G4VITReactionProcess
{
public:
    G4DNAMakeReaction();
    explicit G4DNAMakeReaction(G4VDNAReactionModel*);
    ~G4DNAMakeReaction() override = default;
    G4DNAMakeReaction(const G4DNAMakeReaction& other) = delete;
    G4DNAMakeReaction& operator=(const G4DNAMakeReaction& other) = delete;
    G4bool TestReactibility(const G4Track&,
                            const G4Track&,
                            G4double currentStepTime,
                            G4bool userStepTimeLimit) override;

    std::vector<std::unique_ptr<G4ITReactionChange>> FindReaction(G4ITReactionSet*,
    		const G4double, const G4double, const G4bool) override;

    std::unique_ptr<G4ITReactionChange> MakeReaction(const G4Track&, const G4Track&) override;
    void SetReactionModel(G4VDNAReactionModel*);
    void UpdatePositionForReaction(G4Track&, G4Track&);
    void SetTimeStepComputer(G4VITTimeStepComputer*);
protected:
    const G4DNAMolecularReactionTable*& fMolReactionTable;
    G4VDNAReactionModel* fpReactionModel;
    G4VITTimeStepComputer* fpTimeStepper;
    G4double fTimeStep;
};
#endif