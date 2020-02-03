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
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4VITStepModel.hh"

G4VITStepModel::G4VITStepModel(const G4String& aName)
    : G4VITStepModel(nullptr, nullptr, aName)
{
}

G4VITStepModel::G4VITStepModel(std::unique_ptr<G4VITTimeStepComputer> pTimeStepper,
                               std::unique_ptr<G4VITReactionProcess> pReactionProcess,
                               const G4String& aName)
    : fName(aName)
    , fpTimeStepper(std::move(pTimeStepper))
    , fpReactionProcess(std::move(pReactionProcess))
    , fpReactionTable(nullptr)
    , fType1(-1)
    , fType2(-1)
{
}

void G4VITStepModel::GetApplicable(G4ITType& type1, G4ITType& type2)
{
    type1 = fType1;
    type2 = fType2;
    PrintInfo();
}

void G4VITStepModel::Initialize()
{
    fpReactionProcess->SetReactionTable(fpReactionTable);
    fpTimeStepper->SetReactionTable(fpReactionTable);
    fpTimeStepper->Initialize();
    fpReactionProcess->Initialize();
}

void G4VITStepModel::PrepareNewTimeStep()
{
    fpTimeStepper->Prepare();
}

void G4VITStepModel::SetReactionTable(G4ITReactionTable* pReactionTable)
{
    fpReactionTable = pReactionTable;
}

const G4ITReactionTable* G4VITStepModel::GetReactionTable()
{
    return fpReactionTable ;
}

G4VITTimeStepComputer* G4VITStepModel::GetTimeStepper()
{
    return fpTimeStepper.get();
}

G4VITReactionProcess* G4VITStepModel::GetReactionProcess()
{
    return fpReactionProcess.get();
}

const G4String& G4VITStepModel::GetName()
{
    return fName;
}
