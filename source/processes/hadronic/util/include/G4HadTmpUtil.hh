#ifndef G4HadTmpUtil_hh
#define G4HadTmpUtil_hh 1

#include "globals.hh"
#include <vector>

// Waiting for permission to use namespaces in geant4.....
// namespace G4Had_Tmp_Util
// {
G4int G4lrint(double ad);
G4int G4lint(double ad);
G4int G4rint(double ad);
G4String G4inttostring(int ai);
  template <class A> class G4SortHelperPtr
  {
    public:
    G4SortHelperPtr(A * aA){theA = aA;}
    G4bool operator<(G4SortHelperPtr<A> right) const
    { return *theA < (*right); }
    A & operator * () {return *theA;}
    private:
    A * theA;
  };
#include <algorithm>
template<class A> void G4PtrSort(std::vector<A *> * aList)
{
  std::vector<G4SortHelperPtr<A> > helper;
  for(size_t i=0; i<aList->size(); i++) helper.push_back(aList->operator[](i));
  std::sort(helper.begin(), helper.end());
  for(size_t i=0; i<helper.size(); i++) aList->operator[](i)=&(*helper[i]);
};

// }

#endif
