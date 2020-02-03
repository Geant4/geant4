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
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 


#pragma once

#include <G4VITReactionProcess.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;

/**
  * G4DNAMolecularReaction is the reaction process
  * used in G4DNAMolecularStepByStepModel.
  * It test whether molecules can react after testing.
  * If so, the reaction is made.
  * \deprecated This class will be removed
  */
class G4DNAMolecularReaction : public G4VITReactionProcess
{
public:
    G4DNAMolecularReaction();
    explicit G4DNAMolecularReaction(G4VDNAReactionModel*);
    ~G4DNAMolecularReaction() override = default;
    G4DNAMolecularReaction(const G4DNAMolecularReaction& other) = delete;
    G4DNAMolecularReaction& operator=(const G4DNAMolecularReaction& other) = delete;

    G4bool TestReactibility(const G4Track&,
                            const G4Track&,
                            double currentStepTime,
                            bool userStepTimeLimit) override;

    std::unique_ptr<G4ITReactionChange> MakeReaction(const G4Track&, const G4Track&) override;

    void SetReactionModel(G4VDNAReactionModel*);

protected:
    const G4DNAMolecularReactionTable*& fMolReactionTable;
    G4VDNAReactionModel* fpReactionModel;
};
