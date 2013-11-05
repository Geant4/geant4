#ifndef G4BiasingOperationManager_hh
#define G4BiasingOperationManager_hh 1

class G4VBiasingOperation;
#include <map>
#include <vector>

class G4BiasingOperationManager {
public:
  static G4BiasingOperationManager*               GetInstance();
  const std::vector< G4VBiasingOperation* > GetBiasingOperations()  {return fBiasingOperationVector;}
  G4VBiasingOperation*                       GetBiasingOperation(std::size_t optionID);

public:
  ~G4BiasingOperationManager();
  std::size_t Register(G4VBiasingOperation*);
  

private:
  G4BiasingOperationManager();
  static G4BiasingOperationManager*                                          fInstance;
  static std::vector< G4VBiasingOperation* >                      fBiasingOperationVector;
  static std::map   < G4VBiasingOperation*, std::size_t > fBiasingOperationIDtoPointerMap;
};

#endif
