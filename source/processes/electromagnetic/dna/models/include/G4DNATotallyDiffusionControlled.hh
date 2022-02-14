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
// 20/2/2019
// Author : HoangTRAN

#ifndef G4DNATotallyDiffusionControlled_h
#define G4DNATotallyDiffusionControlled_h 1
#include "G4DNAReactionTypeManager.hh"
#include "G4VReactionType.hh"

class G4MolecularConfiguration;

class G4DNAMolecularReactionTable;
class G4VDNAReactionModel;

class G4DNATotallyDiffusionControlled 
    : public G4VReactionType
{
public:
    G4DNATotallyDiffusionControlled();
    ~G4DNATotallyDiffusionControlled() override;

    G4DNATotallyDiffusionControlled(const G4DNATotallyDiffusionControlled& other) = delete;
    G4DNATotallyDiffusionControlled& operator=(const G4DNATotallyDiffusionControlled& other) = delete;

    G4double GetTimeToEncounter(const G4Track& trackA,
                              const G4Track& trackB) override;

    G4bool GeminateRecombinationProbability(const G4MolecularConfiguration* pConfMolA,
                                          const G4MolecularConfiguration* pConfMolB) override;
private:
    G4double GetDiffusionCoefficient(const G4MolecularConfiguration* pMA,
                                   const G4MolecularConfiguration* pMB);
};
#endif