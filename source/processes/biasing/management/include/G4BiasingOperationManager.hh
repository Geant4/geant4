#ifndef G4BiasingOperationManager_hh
#define G4BiasingOperationManager_hh 1

class G4VBiasingOperation;
#include <map>
#include <vector>
#include "G4Cache.hh"
#include "G4ThreadLocalSingleton.hh"

//This class is a thread-local singleton 
class G4BiasingOperationManager {
    friend class G4ThreadLocalSingleton<G4BiasingOperationManager>;
public:
  static G4BiasingOperationManager*                  GetInstance();
  const std::vector< G4VBiasingOperation* > GetBiasingOperations() {return fBiasingOperationVector.Get();}
  G4VBiasingOperation*                       GetBiasingOperation(std::size_t optionID);

public:
  ~G4BiasingOperationManager();
  std::size_t Register(G4VBiasingOperation*);
  
  
private:
  G4BiasingOperationManager();
  static G4VectorCache<G4VBiasingOperation*> fBiasingOperationVector;
  static G4MapCache<G4VBiasingOperation*,std::size_t > fBiasingOperationIDtoPointerMap;
  //static G4BiasingOperationManager*                                          fInstance;
  //static std::vector< G4VBiasingOperation* >                      fBiasingOperationVector;
  //static std::map   < G4VBiasingOperation*, std::size_t > fBiasingOperationIDtoPointerMap;
};

#endif
