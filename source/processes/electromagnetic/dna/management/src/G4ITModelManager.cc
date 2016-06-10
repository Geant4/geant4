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
// $Id: G4ITModelManager.cc 87375 2014-12-02 08:17:28Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr)
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4ITModelManager.hh"
#include "G4ITType.hh"
#include "G4UnitsTable.hh"

using namespace std;

G4ITModelManager::G4ITModelManager() :
    fIsInitialized(FALSE)
{
  ;
}

G4ITModelManager::~G4ITModelManager()
{
  //dtor
  mapModels::iterator it;

  for (it = fModels.begin(); it != fModels.end(); it++)
  {
    delete it->second;
  }
  fModels.clear();
}

G4ITModelManager::G4ITModelManager(const G4ITModelManager& right)
{
  mapModels::const_iterator it = right.fModels.begin();

  for (; it != right.fModels.end(); it++)
  {
    fModels[it->first] = it->second->Clone();
  }

  fIsInitialized = right.fIsInitialized;
}

G4ITModelManager& G4ITModelManager::operator=(const G4ITModelManager& rhs)
{
  if (this == &rhs) return *this; // handle self assignment
  //assignment operator
  return *this;
}

void G4ITModelManager::Initialize()
{
  mapModels::iterator it = fModels.begin();

  for (; it != fModels.end(); it++)
  {
    G4VITStepModel* model = it->second;
    if (model != 0)
    {
      model->Initialize();
    }
  }
}

void G4ITModelManager::SetModel(G4VITStepModel* aModel, G4double startingTime)
{
  assert(fIsInitialized == FALSE);
  if (fIsInitialized == true)
  {
    G4ExceptionDescription exceptionDescription;
    exceptionDescription
        << "You are trying to insert a new model after initialization of th model manager.";
    G4Exception("G4ITModelManager::SetModel", "ITModelManager001",
                FatalErrorInArgument, exceptionDescription);
  }
  fModels[startingTime] = aModel;
}

G4VITStepModel* G4ITModelManager::GetModel(const G4double globalTime)
{
  if (!fModels.empty())
  {
    mapModels::reverse_iterator rit = fModels.rbegin();
    if (rit != fModels.rend())
    {
      if (globalTime > rit->first)
      {
        return rit->second;
      }
      else
      {
        mapModels::iterator it = fModels.begin();

        if (globalTime < it->first)
        {
          G4ExceptionDescription exceptionDescription;
          exceptionDescription << "No model was found at time ";
          exceptionDescription << G4BestUnit(globalTime, "Time");
          exceptionDescription << ". The first model is registered at time : ";
          exceptionDescription << G4BestUnit(it->first, "Time") << ". ";
          G4Exception("G4ITModelManager::GetModel", "ITModelManager003",
                      FatalErrorInArgument, exceptionDescription);
        }

        it = fModels.lower_bound(globalTime);

        if (it != fModels.end()) return it->second;
      }
    }
  }

  G4ExceptionDescription exceptionDescription;
  exceptionDescription << "No model was found.";
  G4Exception("G4ITModelManager::GetModel", "ITModelManager004",
              FatalErrorInArgument, exceptionDescription);
  return 0;
}
