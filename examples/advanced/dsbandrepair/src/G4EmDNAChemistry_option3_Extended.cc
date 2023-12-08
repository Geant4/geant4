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
/// \file G4EmDNAChemistry_option3_Extended.cc
/// \brief Implementation of the G4EmDNAChemistry_option3_Extended class

#include "G4EmDNAChemistry_option3_Extended.hh"

#include "G4DNAMolecule.hh"
#include "G4MoleculeTable.hh"

// particles
#include "G4H2O.hh"
#include "G4H2.hh"
#include "G4Hydrogen.hh"
#include "G4OH.hh"
#include "G4H3O.hh"
#include "G4Electron_aq.hh"

#include "G4H2O2.hh"
#include "G4O2.hh"
#include "G4HO2.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAIndependentReactionTimeModel.hh"
#include "G4DNAMolecularReactionTable.hh"

#include "G4ChemicalMoleculeFinder.hh"
// factory
#include "G4PhysicsConstructorFactory.hh"

G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAChemistry_option3_Extended);

void G4EmDNAChemistry_option3_Extended::ConstructParticle()
{
    ConstructMolecule(); //from G4EmDNAChemistry_option3
    //dna molecules 
    G4Deoxyribose::Definition();
    G4Phosphate::Definition();
    G4Adenine::Definition();
    G4Guanine::Definition();
    G4Thymine::Definition();
    G4Cytosine::Definition();
    G4Histone::Definition();
    //damaged molecules
    G4DamagedDeoxyribose::Definition();
    G4DamagedAdenine::Definition();
    G4DamagedGuanine::Definition();
    G4DamagedThymine::Definition();
    G4DamagedCytosine::Definition();
    G4ModifiedHistone::Definition();
    //________________DNA_______________________________________________
    auto table = G4MoleculeTable::Instance();
    table->CreateConfiguration("Deoxyribose",G4Deoxyribose::Definition());
    table->CreateConfiguration("Phosphate",G4Phosphate::Definition());
    table->CreateConfiguration("Adenine",G4Adenine::Definition());
    table->CreateConfiguration("Thymine",G4Thymine::Definition());
    table->CreateConfiguration("Guanine",G4Guanine::Definition());
    table->CreateConfiguration("Cytosine",G4Cytosine::Definition());
    table->CreateConfiguration("Histone",G4Histone::Definition());
  
    //damaged DNAElement Configuration

    table->CreateConfiguration("Damaged_Deoxyribose",
    G4DamagedDeoxyribose::Definition());
    table->CreateConfiguration("Damaged_Adenine",
    G4DamagedAdenine::Definition());
    table->CreateConfiguration("Damaged_Thymine",
    G4DamagedThymine::Definition());
    table->CreateConfiguration("Damaged_Guanine",
    G4DamagedGuanine::Definition());
    table->CreateConfiguration("Damaged_Cytosine",
    G4DamagedCytosine::Definition());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAChemistry_option3_Extended::ConstructReactionTable(
    G4DNAMolecularReactionTable* theReactionTable)
{
    G4EmDNAChemistry_option3::ConstructReactionTable(theReactionTable);//from G4EmDNAChemistry_option3

    //Get the molecular configuration
    auto table = G4MoleculeTable::Instance();
    G4MolecularConfiguration* OH = table->GetConfiguration("OH");
    G4MolecularConfiguration* OHm = table->GetConfiguration("OHm");
    G4MolecularConfiguration* e_aq = table->GetConfiguration("e_aq");
    G4MolecularConfiguration* H2 = table->GetConfiguration("H2");
    G4MolecularConfiguration* H3Op = table->GetConfiguration("H3Op");
    G4MolecularConfiguration* H = table->GetConfiguration("H");
    G4MolecularConfiguration* H2O2 = table->GetConfiguration("H2O2");

    // DNA additions--------------------------------------------------
    G4MolecularConfiguration* deoxyribose = table->GetConfiguration("Deoxyribose");
    G4MolecularConfiguration* adenine = table->GetConfiguration("Adenine");
    G4MolecularConfiguration* guanine = table->GetConfiguration("Guanine");
    G4MolecularConfiguration* thymine = table->GetConfiguration("Thymine");
    G4MolecularConfiguration* cytosine = table->GetConfiguration("Cytosine");
    G4MolecularConfiguration* histone = table->GetConfiguration("Histone");

    G4MolecularConfiguration* damage_deoxyribose = table->GetConfiguration("Damaged_Deoxyribose");
    G4MolecularConfiguration* damage_adenine = table->GetConfiguration("Damaged_Adenine");
    G4MolecularConfiguration* damage_guanine = table->GetConfiguration("Damaged_Guanine");
    G4MolecularConfiguration* damage_thymine = table->GetConfiguration("Damaged_Thymine");
    G4MolecularConfiguration* damage_cytosine = table->GetConfiguration("Damaged_Cytosine");


    // OH and DNA

    // 2-Deoxyribose + OH -> damagedDeoxyribose
    G4DNAMolecularReactionData* reactionData = new G4DNAMolecularReactionData(
    1.80e9*(1e-3*m3/(mole*s)), deoxyribose, OH);
    reactionData->AddProduct(damage_deoxyribose);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // adenine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    6.10e9*(1e-3*m3/(mole*s)), adenine, OH);
    reactionData->AddProduct(damage_adenine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // guanine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    9.20e9*(1e-3*m3/(mole*s)), guanine, OH);
    reactionData->AddProduct(damage_guanine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // thymine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    6.40e9*(1e-3*m3/(mole*s)), thymine, OH);
    reactionData->AddProduct(damage_thymine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // cytosine + OH -> ...
    reactionData = new G4DNAMolecularReactionData(
    6.10e9*(1e-3*m3/(mole*s)), cytosine, OH);
    reactionData->AddProduct(damage_cytosine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // Hydrated e- and DNA

    // Deoxyribose + Hydrated e-  -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.01e9*(1e-3*m3/(mole*s)), deoxyribose, e_aq);
    reactionData->AddProduct(damage_deoxyribose);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // adenine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    9e9*(1e-3*m3/(mole*s)), adenine, e_aq);
    reactionData->AddProduct(damage_adenine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // guanine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    14e9*(1e-3*m3/(mole*s)), guanine, e_aq);
    reactionData->AddProduct(damage_guanine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // thymine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    18e9*(1e-3*m3/(mole*s)), thymine, e_aq);
    reactionData->AddProduct(damage_thymine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // cytosine + Hydrated e- -> ...
    reactionData = new G4DNAMolecularReactionData(
    13e9*(1e-3*m3/(mole*s)), cytosine, e_aq);
    reactionData->AddProduct(damage_cytosine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // Radical H and DNA

    // Deoxyribose + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.029e9*(1e-3*m3/(mole*s)), deoxyribose, H);
    reactionData->AddProduct(damage_deoxyribose);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // adenine + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.10e9*(1e-3*m3/(mole*s)), adenine, H);
    reactionData->AddProduct(damage_adenine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // thymine + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.57e9*(1e-3*m3/(mole*s)), thymine, H);
    reactionData->AddProduct(damage_thymine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    // cytosine + Radical H -> ...
    reactionData = new G4DNAMolecularReactionData(
    0.092e9*(1e-3*m3/(mole*s)), cytosine, H);
    reactionData->AddProduct(damage_cytosine);
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    //histone + all molecules -> modification(or "damage")

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, OH);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4OH::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, OHm);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4OH::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, e_aq);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4Electron_aq::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H2);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4H2::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H3Op);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4H3O::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4Hydrogen::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);

    reactionData = new G4DNAMolecularReactionData(
    0.0*(1e-3*m3/(mole*s)), histone, H2O2);
    reactionData->AddProduct(histone);
    reactionData->SetEffectiveReactionRadius(
    2.4*nm + G4H2O2::Definition()->GetVanDerVaalsRadius());
    reactionData->SetReactionType(1);
    theReactionTable->SetReaction(reactionData);
}

