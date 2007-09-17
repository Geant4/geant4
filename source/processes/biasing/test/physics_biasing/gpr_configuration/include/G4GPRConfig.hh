#ifndef G4GPRCONFIG_HH
#define G4GPRCONFIG_HH

#include <vector>
#include <map>

template <typename ProcessSelection, typename ParticleSelection> 
struct G4GPRConfig : public ProcessSelection, public ParticleSelection {

  public:
    
  typedef std::vector<unsigned> ProcessList;
  typedef std::vector<G4ParticleDefinition*> ParticleList;
  typedef std::map<G4ParticleDefinition*, ProcessList> Result;

  template <typename List>
  Result GetConfig() const
  {
    Result result;

    ParticleList particleList = ParticleSelection::GetList();

    ParticleList::iterator iter = particleList.begin();

    while (iter != particleList.end()) {
      ProcessList processList = ProcessSelection::template GetList<List>(*iter);
      result[*iter] = processList;

      G4cout<<"jane processlist for "<<(*iter)->GetParticleName();
      ProcessList::iterator iter2= processList.begin();
      
      while(iter2 != processList.end()) {
	G4cout<<*iter2<<" ";
	iter2++;
      }
      
      G4cout<<G4endl;
      iter++;
    }

    return result;

  }

};

#endif
