#include "G4VBiasingOperation.hh"
#include "G4BiasingOperationManager.hh"

G4VBiasingOperation::G4VBiasingOperation(G4String name)
  : fName(name)
{
  // -- have to move from initialization to here, as giving error on windows:
  fUniqueID = G4BiasingOperationManager::GetInstance()->Register(this);
}

G4VBiasingOperation::~G4VBiasingOperation()
{}
