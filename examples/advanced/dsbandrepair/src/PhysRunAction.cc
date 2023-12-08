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
/// \file PhysRunAction.cc
/// \brief Implementation of the PhysRunAction class

#include "PhysRunAction.hh"
#include "G4Run.hh"
#include "PhysAnalysis.hh"
#include "InformationKeeper.hh"
#include "G4DNAChemistryManager.hh"
#ifdef USE_MPI
#include "G4MPImanager.hh"
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

PhysRunAction::PhysRunAction() : G4UserRunAction()
{
    CreateNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysRunAction::BeginOfRunAction(const G4Run*)
{
    PhysAnalysis* myana = PhysAnalysis::GetAnalysis();
    std::ostringstream fname;
    G4String slash = "";
#if defined(_WIN32) || defined(WIN32)
    slash= "\\";
#else
    G4String slashu = "/";
    slash = slashu;
#endif

#ifdef USE_MPI
    G4int rank = G4MPImanager::GetManager()->GetRank();
    fname<<InformationKeeper::Instance()->GetPhysOutFolderName()<<slash<<"phys_output"<<"_rank"<<rank;
#else
    fname<<InformationKeeper::Instance()->GetPhysOutFolderName()<<slash<<"phys_output";
#endif
    myana->OpenFile(fname.str());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysRunAction::EndOfRunAction(const G4Run*)
{   
    WriteNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysRunAction::CreateNtuple()
{
    auto myana = PhysAnalysis::GetAnalysis();
    myana->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void PhysRunAction::WriteNtuple()
{
    PhysAnalysis* myana = PhysAnalysis::GetAnalysis();
    myana-> Save();
	myana-> Close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






