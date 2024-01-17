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

#include "G4ITModelManager.hh"
#include "G4VITStepModel.hh"
#include "G4ITType.hh"
#include "G4UnitsTable.hh"


G4ITModelManager::G4ITModelManager()
     
= default;

G4ITModelManager::~G4ITModelManager() = default;

void G4ITModelManager::Initialize()
{
    std::sort(fModelInfoList.begin(), fModelInfoList.end(),
              [](const ModelInfo& lhs, const ModelInfo& rhs) {
                  return lhs.fStartingTime < rhs.fStartingTime;
              });

    for (const auto& modelInfo : fModelInfoList)
    {
        modelInfo.fpModel->Initialize();
    }

    fIsInitialized = true;
}

void G4ITModelManager::SetModel(G4VITStepModel* pModel,
                                const G4double startingTime,
                                const G4double endTime)
{
    assert(pModel != nullptr);
    if (fIsInitialized)
    {
        G4ExceptionDescription exceptionDescription;
        exceptionDescription
            << "You are trying to insert a new model after initializing the model manager.";
        G4Exception("G4ITModelManager::SetModel", "ITModelManager001",
                    FatalErrorInArgument, exceptionDescription);
    }

    fModelInfoList.emplace_back(ModelInfo({ startingTime, endTime, std::unique_ptr<G4VITStepModel>(pModel) }));
}

std::vector<G4VITStepModel*> G4ITModelManager::GetActiveModels(G4double globalTime) const
{
    std::vector<G4VITStepModel*> activeModels;

    for (const auto& modelInfo : fModelInfoList)
    {
        if (modelInfo.fStartingTime < globalTime && modelInfo.fEndTime > globalTime)
        {
            activeModels.push_back(modelInfo.fpModel.get());
        }
    }

    return activeModels;
}

