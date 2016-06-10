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
// $Id: G4ITModelHandler.cc 85244 2014-10-27 08:24:13Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITModelHandler.hh"
#include <assert.h>

G4ITModelHandler::G4ITModelHandler()
{
  //ctor
  fIsInitialized = false;
  fTimeStepComputerFlag = false;
  fReactionProcessFlag = false;

  size_t IT_size(G4ITType::size());

  fModelManager.assign(IT_size, std::vector<G4ITModelManager*>());
  for (G4int i = 0; i < (int) IT_size; i++)
  {
    fModelManager[i].assign(IT_size, 0);
  }
}

G4ITModelHandler::~G4ITModelHandler()
{
  //dtor
  G4int size = fModelManager.size();

  for (G4int i = 0; i < size; i++)
  {
    for (G4int j = 0; j <= i; j++)
    {
      if (fModelManager[i][j])
      {
        delete fModelManager[i][j];
        fModelManager[i][j] = 0;
        fModelManager[j][i] = 0;
      }
    }
  }
  fModelManager.clear();
}

G4ITModelHandler::G4ITModelHandler(const G4ITModelHandler& other)
{
  //copy ctor
  size_t IT_size(G4ITType::size());

  fModelManager.assign(IT_size, std::vector<G4ITModelManager*>());
  for (int i = 0; i < (int) IT_size; i++)
  {
    fModelManager[i].assign(IT_size, 0);
    for (int j = 0; j < (int) IT_size; j++)
    {
      if (other.fModelManager[i][j] != 0)
      {
        fModelManager[i][j] = new G4ITModelManager(
            *(other.fModelManager[i][j]));
      }
    }
  }

  fIsInitialized = other.fIsInitialized;
  fTimeStepComputerFlag = other.fTimeStepComputerFlag;
  fReactionProcessFlag = other.fReactionProcessFlag;
}

G4ITModelHandler& G4ITModelHandler::operator=(const G4ITModelHandler& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

void G4ITModelHandler::Initialize()
{
  fIsInitialized = true;

  for (G4int i = 0; i < int(fModelManager.size()); i++)
  {
    for (G4int j = 0; j <= i; j++)
    {
      G4ITModelManager* modman = fModelManager[i][j];
      if (modman)
      {
        modman->Initialize();
      }
    }
  }
}

void G4ITModelHandler::RegisterModel(G4VITStepModel* aModel,
                                     G4double startingTime)
{
  assert(aModel != 0);

  //________________________________________________
  // Prepare the correct model manager
  if (fModelManager.empty())
  {
    size_t IT_size(G4ITType::size());
    fModelManager.assign(IT_size, std::vector<G4ITModelManager*>());

    for (int i = 0; i < (int) IT_size; i++)
    {
      fModelManager[i].assign((size_t) i, 0);
    }
  }

  G4ITType type1;
  G4ITType type2;

  //________________________________________________
  // Retrieve applicability
  aModel->IsApplicable(type1, type2);

  if (type1 > type2)
  {
    G4ITType buffer(-1);
    buffer = type1;
    type1 = type2;
    type2 = buffer;
  }

  if (fModelManager[type1][type2] == 0)
  {
    fModelManager[type1][type2] = new G4ITModelManager();
  }

  fModelManager[type1][type2]->SetModel(aModel, startingTime);

  //________________________________________________
  // Setup ITStepManager
  if (aModel->GetTimeStepper())
  {
    fTimeStepComputerFlag = true;
  }
  if (aModel->GetReactionProcess())
  {
    fReactionProcessFlag = true;
  }
}

void G4ITModelHandler::SetModel(G4ITType type1,
                                G4ITType type2,
                                G4VITStepModel* aModel,
                                G4double startingTime)
{
  assert(aModel == 0);

  if (type1 > type2)
  {
    G4ITType buffer(-1);
    buffer = type1;
    type1 = type2;
    type2 = buffer;
  }

  if (type1 > (int) fModelManager.capacity())
  {
    fModelManager.reserve(type1);
  }

  if (type2 > (int) fModelManager[type1].capacity())
  {
    fModelManager[type1].reserve(type2);
  }

  fModelManager[type1][type2]->SetModel(aModel, startingTime);
}

G4VITStepModel* G4ITModelHandler::GetModel(G4ITType type1,
                                           G4ITType type2,
                                           const G4double globalTime)
{
  if (fModelManager.empty())
  {
    return 0;
  }

  if ((int) fModelManager.size() < type1) return 0;

  std::vector<G4ITModelManager*>* v = &(fModelManager.at(type1));

  if ((int) v->size() < type2) return 0;

  if (v->at(type2))
  {
    return v->at(type2)->GetModel(globalTime);
  }
  return 0;
}
