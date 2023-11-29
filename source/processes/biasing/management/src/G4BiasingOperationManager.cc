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
#include "G4BiasingOperationManager.hh"

//G4BiasingOperationManager*                       G4BiasingOperationManager::fInstance = 0;
G4VectorCache< G4VBiasingOperation* >              G4BiasingOperationManager::fBiasingOperationVector;
G4MapCache< G4VBiasingOperation*, std::size_t > G4BiasingOperationManager::fBiasingOperationIDtoPointerMap;

G4BiasingOperationManager::G4BiasingOperationManager()
{}

G4BiasingOperationManager::~G4BiasingOperationManager()
{}

G4BiasingOperationManager* G4BiasingOperationManager::GetInstance()
{
    //Create an instance for each thread.
    static G4ThreadLocalSingleton<G4BiasingOperationManager> instance;
    return instance.Instance();
//  if (fInstance == 0) fInstance = new G4BiasingOperationManager();
//  return fInstance;
}

std::size_t G4BiasingOperationManager::Register(G4VBiasingOperation* option)
{
  std::size_t optionUniqueID = fBiasingOperationVector.Size();
  
  fBiasingOperationVector.Push_back(option);
  fBiasingOperationIDtoPointerMap[option] = optionUniqueID;
  
  return optionUniqueID;
}

G4VBiasingOperation* G4BiasingOperationManager::GetBiasingOperation(std::size_t optionID)
{
  if (optionID < fBiasingOperationVector.Size())
    return fBiasingOperationVector[(G4int)optionID];
  else
    return nullptr;
}
