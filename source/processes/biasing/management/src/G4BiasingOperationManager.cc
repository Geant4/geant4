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
  if (optionID < fBiasingOperationVector.Size()) return fBiasingOperationVector[optionID];
  else return 0;
}




