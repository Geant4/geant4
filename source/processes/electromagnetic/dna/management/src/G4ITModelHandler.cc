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

#include "G4ITModelHandler.hh"
#include "G4ITModelManager.hh"
#include "G4VITStepModel.hh"
#include <assert.h>

G4ITModelHandler::G4ITModelHandler()
{
    fIsInitialized = false;
    fTimeStepComputerFlag = false;
    fReactionProcessFlag = false;
}

G4ITModelHandler::~G4ITModelHandler() = default;

void G4ITModelHandler::Initialize()
{
    fpModelManager->Initialize();
    fIsInitialized = true;
}

void G4ITModelHandler::RegisterModel(G4VITStepModel* pModel,
                                     G4double startingTime)
{
    if(fFinalize)
    {
      return;
    }
    assert(pModel != nullptr);

    G4ITType type1;
    G4ITType type2;

    pModel->GetApplicable(type1, type2);

    if (type1 != type2)
    {
        G4Exception("G4ITModelHandler::RegisterModel",
                    "FeatureDisabled",
                    FatalException,
                    "Models for different type ids are not supported anymore. This feature will be superseded.");
    }

    if(!fpModelManager)
    {
        fpModelManager.reset(new G4ITModelManager());
    }

    fpModelManager->SetModel(pModel, startingTime);

    //________________________________________________
    // Setup ITStepManager
    if (pModel->GetTimeStepper())
    {
        fTimeStepComputerFlag = true;
    }
    if (pModel->GetReactionProcess())
    {
        fReactionProcessFlag = true;
    }
}

std::vector<G4VITStepModel*> G4ITModelHandler::GetActiveModels(G4double globalTime) const
{
    if(!fpModelManager)
    {
        return {};
    }
    return fpModelManager->GetActiveModels(globalTime);
}

bool G4ITModelHandler::GetTimeStepComputerFlag()
{
    return fTimeStepComputerFlag;
}

bool G4ITModelHandler::GetReactionProcessFlag()
{
    return fReactionProcessFlag;
}

void G4ITModelHandler::Reset()
{
    fpModelManager.reset();
}

void G4ITModelHandler::Finalize()
{
  fFinalize = true;
}