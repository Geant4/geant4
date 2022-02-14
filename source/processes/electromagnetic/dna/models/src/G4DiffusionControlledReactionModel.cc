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
// Author: Hoang TRAN

#include "G4DiffusionControlledReactionModel.hh"
#include "G4Track.hh"
#include "G4DNAMolecularReactionTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Exp.hh"
#include "G4IRTUtils.hh"
#include "G4SystemOfUnits.hh"
#include "G4DNAReactionTypeManager.hh"
#include "G4Electron_aq.hh"
G4DiffusionControlledReactionModel::G4DiffusionControlledReactionModel()
    : G4VDNAReactionModel()
    , fpReactionData(nullptr)
    , fReactionTypeManager(nullptr)
{
}

G4DiffusionControlledReactionModel::~G4DiffusionControlledReactionModel() = default;

void G4DiffusionControlledReactionModel::Initialise(const G4MolecularConfiguration* pMolecule,
                                                const G4Track&)
{
    fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

void G4DiffusionControlledReactionModel::InitialiseToPrint(const G4MolecularConfiguration* pMolecule)
{
    fpReactionData = fpReactionTable->GetReactionData(pMolecule);
}

G4double G4DiffusionControlledReactionModel::GetReactionRadius(const G4MolecularConfiguration* pMol1,
                                                               const G4MolecularConfiguration* pMol2)
{
    auto reactionData = fpReactionTable->GetReactionData(pMol1, pMol2);
    if(reactionData == nullptr)
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription << "No reactionData"
                             <<" for : "<<pMol1->GetName()
                             <<" and "<<pMol2->GetName();
        G4Exception("G4DiffusionControlledReactionModel"
                    "::GetReactionRadius()", "G4DiffusionControlledReactionModel00",
                    FatalException, exceptionDescription);
    }
    G4double kobs = reactionData->GetObservedReactionRateConstant();
    G4double D;
    if(pMol1 == pMol2)
    {    
        D = (pMol1->GetDiffusionCoefficient());
    }
    else
    {
        D = (pMol1->GetDiffusionCoefficient() +
             pMol2->GetDiffusionCoefficient());//
    }
    
    if ( D == 0 )
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription << "D = "<< D
        << " is uncorrected"
        <<" for : "<<pMol1->GetName()
        <<" and "<<pMol2->GetName();
        G4Exception("G4DiffusionControlledReactionModel"
                    "::GetReactionRadius()", "G4DiffusionControlledReactionModel01",
                    FatalException, exceptionDescription);
    }

    G4double Reff = kobs / ( 4 * CLHEP::pi * D * Avogadro );
    return Reff;
}

G4double G4DiffusionControlledReactionModel::GetReactionRadius(G4int i)
{
    auto pMol1 = (*fpReactionData)[i]->GetReactant1();
    auto pMol2 = (*fpReactionData)[i]->GetReactant2();
    return GetReactionRadius(pMol1,pMol2);
}

void G4DiffusionControlledReactionModel::SetReactionTypeManager(G4VReactionTypeManager* typeManager)
{
    fReactionTypeManager = ((G4DNAReactionTypeManager*)typeManager);
}
