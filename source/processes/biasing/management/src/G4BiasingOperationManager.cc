#include "G4BiasingOperationManager.hh"

G4BiasingOperationManager*                       G4BiasingOperationManager::fInstance = 0;
std::vector< G4VBiasingOperation* >              G4BiasingOperationManager::fBiasingOperationVector;
std::map   < G4VBiasingOperation*, std::size_t > G4BiasingOperationManager::fBiasingOperationIDtoPointerMap;

G4BiasingOperationManager::G4BiasingOperationManager()
{}

G4BiasingOperationManager::~G4BiasingOperationManager()
{}

G4BiasingOperationManager* G4BiasingOperationManager::GetInstance()
{
  if (fInstance == 0) fInstance = new G4BiasingOperationManager();
  return fInstance;
}

std::size_t G4BiasingOperationManager::Register(G4VBiasingOperation* option)
{
  std::size_t optionUniqueID = fBiasingOperationVector.size();
  
  fBiasingOperationVector.push_back(option);
  fBiasingOperationIDtoPointerMap[option] = optionUniqueID;
  
  return optionUniqueID;
}

G4VBiasingOperation* G4BiasingOperationManager::GetBiasingOperation(std::size_t optionID)
{
  if (optionID < fBiasingOperationVector.size()) return fBiasingOperationVector[optionID];
  else return 0;
}


