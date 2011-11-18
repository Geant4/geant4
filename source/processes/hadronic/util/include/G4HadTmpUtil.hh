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
#ifndef G4HadTmpUtil_hh
#define G4HadTmpUtil_hh 1

#include "globals.hh"
#include <vector>

// Waiting for permission to use namespaces in geant4.....
// namespace G4Had_Tmp_Util
// {
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
}

struct G4Delete { template<class T> void operator() (T * aT) {delete aT;} };

// }

#endif
