#ifndef G4GPRPHYSICSLISTMANAGER_HH
#define G4GPRPHYSICSLISTMANAGER_HH

#include "G4GPRPhysicsList.hh"
#include <vector>

class G4GPRPhysicsListManager {

public:
  
  G4GPRPhysicsListManager()
    :pDefaultList(new G4GPRPhysicsList("Default"))
    ,pActiveList(0)
    ,fListChanged(false)
    ,fNActiveLists(0)
  {
    fLists.push_back(pDefaultList);
  }
  
  G4GPRPhysicsList* GetPhysicsList(const G4String& name)
  {
    // jane fixme
    if (name == "") return pDefaultList;
    std::vector<G4GPRPhysicsList*>::iterator iter = fLists.begin();
    
    while (iter != fLists.end()) {
      if ((*iter)->GetName() == name) return *iter;
      iter++;
    }
    
    return 0;
  }
  
  void SetDefaultList(G4GPRPhysicsList* list)
  {
    pDefaultList = list;
    fLists.push_back(list);
    fListChanged = true;
  }

  G4bool ListChanged() {return fListChanged;}

  void ChangeState(G4GPRPhysicsList* list) 
  {    
    fListChanged = true;
    G4cout<<"jane physics list chagned "<<list->GetName()<<" "<<list->GetState()<<G4endl;
    if (list->GetState()) {      
      pActiveList = list; 
      fNActiveLists++;
	
      return;
    }
    
    G4cout<<"jane change state "<<pActiveList<<" "<<list<<G4endl;
    if (pActiveList == list) pActiveList = 0;
    
    return;
      
  }
  void ResetFlags()
  {
    fNActiveLists = 1;
    fListChanged = false;
  }
  
  G4GPRPhysicsList* GetDefaultList() {return pDefaultList;}

  G4String GetActiveListName() 
  {
    return (0 != pActiveList ? pActiveList->GetName() : pDefaultList->GetName());
  }
  
  G4GPRPhysicsList* GetActiveList()
  {
    return (0 != pActiveList ? pActiveList : pDefaultList);
  }
  
  void Register(G4GPRPhysicsList* list) 
  {
    fLists.push_back(list);
    
    list->AddObserver(this, &G4GPRPhysicsListManager::ChangeState);
  }
  
  
private:
  
  G4GPRPhysicsList* pDefaultList;
  G4GPRPhysicsList* pActiveList;
  
  G4bool fListChanged;
  unsigned fNActiveLists;
  
  std::vector<G4GPRPhysicsList*> fLists;
  
};

#endif
