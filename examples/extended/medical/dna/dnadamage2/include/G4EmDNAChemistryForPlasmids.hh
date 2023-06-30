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
/// \file G4EmDNAChemistryForPlasmids.hh
/// \brief Definition of the Chemistry parameters with DNA reactions
/*
 * G4EmDNAChemistryForPlasmids.hh
 *
 *  Created on: Feb 10, 2021
 *      Author: J. Naoki D. Kondo
 *              W. G. Shin, J. Ramos-Mendez and B. Faddegon
*/

#ifndef DNADAMAGE2_ChemistryForPlasmids_hh
#define DNADAMAGE2_ChemistryForPlasmids_hh 1

#include "G4UIcmdWithADouble.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4VUserChemistryList.hh"
#include <G4UImessenger.hh>
#include "globals.hh"

class G4DNAMolecularReactionTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4EmDNAChemistryForPlasmids:public G4VUserChemistryList,
                                  public G4VPhysicsConstructor
{
public:
  G4EmDNAChemistryForPlasmids();
  G4EmDNAChemistryForPlasmids(G4double, G4double);
  ~G4EmDNAChemistryForPlasmids() override = default;

  void ConstructParticle() override
  {
    ConstructMolecule();
  }

  void ConstructMolecule() override;
  void ConstructProcess() override;

  void ConstructDissociationChannels() override;
  void ConstructReactionTable(G4DNAMolecularReactionTable*) override;
  void ConstructTimeStepModel(G4DNAMolecularReactionTable*) override;

  void DeclareDMSOAndDNAReactions(G4DNAMolecularReactionTable*);
  void DeclareOxygenReactions(G4DNAMolecularReactionTable*);

private:
  G4double fDMSO   = 1E-4;
  G4double fOxygen = 0.25E-3;

  //G4DNAMolecularReactionTable* fReactionTable = nullptr;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
