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
 * G4DNAMolecularIRTModel.hh
 *
 *  Created on: Jul 23, 2019
 *      Author: W. G. Shin
 *              J. Ramos-Mendez and B. Faddegon
*/

#pragma once

#include <G4String.hh>
#include <G4VITStepModel.hh>

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;

class G4DNAMolecularIRTModel : public G4VITStepModel
{
public:
    G4DNAMolecularIRTModel(const G4String& name = "DNAMolecularIRTModel");
    G4DNAMolecularIRTModel(const G4String& name,
                                  std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
                                  std::unique_ptr<G4VITReactionProcess> pReactionProcess);
    G4DNAMolecularIRTModel& operator=(const G4DNAMolecularIRTModel&) = delete;
    G4DNAMolecularIRTModel(const G4DNAMolecularIRTModel&) = delete;
    ~G4DNAMolecularIRTModel() override;

    void PrintInfo() override;
    void Initialize() override;

    void SetReactionModel(G4VDNAReactionModel*);
    G4VDNAReactionModel* GetReactionModel();

protected:
    const G4DNAMolecularReactionTable*& fMolecularReactionTable;
    std::unique_ptr<G4VDNAReactionModel> fpReactionModel;
};

