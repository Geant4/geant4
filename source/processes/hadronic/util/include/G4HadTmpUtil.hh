//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
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
  for(size_t j=0; j<helper.size(); j++) aList->operator[](j)=&(*helper[j]);
};

struct G4Delete { template<class T> void operator() (T * aT) {delete aT;} };

// }

#endif
