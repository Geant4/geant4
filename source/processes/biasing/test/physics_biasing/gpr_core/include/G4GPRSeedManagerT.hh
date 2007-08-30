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
// $Id: G4GPRSeedManagerT.hh,v 1.3 2007-08-30 19:37:45 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRSEEDMANAGERT_HH
#define G4GPRSEEDMANAGERT_HH

#include "G4GPRPlacement.hh"
#include "G4GPRManagerT.hh"
#include "G4GPRSeedT.hh"
#include <vector>

template <>
template <typename List>
class G4GPRManagerT< G4GPRSeedT<List> > {

public:

  typedef typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper Wrapper;
  typedef std::vector< G4GPRSeedT<List>* > Store;
  typedef std::vector<typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper> Result;

  void Register(G4GPRSeedT<List>* component)
  {
    fStore.push_back(component);
  }
  
  void operator()(Result*& result) 
  {
    typename Store::iterator iter = fStore.begin();

    std::vector<Wrapper> first;
    std::vector<Wrapper> last;
    std::vector<Wrapper> append;
    std::map<G4int, Wrapper> placed;

    while (iter != fStore.end()) {
      
      if ((*iter)->IsActive()) {
	G4int idx = (*iter)->Placement();
	
	switch (idx) {
	  
	case G4GPRPlacement::First :
	  first.push_back((*iter)->GetWrapper());
	  break;
	case G4GPRPlacement::Last :
	  last.push_back((*iter)->GetWrapper());
	  break;
	case G4GPRPlacement::Append :
	  append.push_back((*iter)->GetWrapper());
	  break;
	default :
	  placed[idx] = (*iter)->GetWrapper();
	}

      }
      iter++;
    }

    assert (first.size() <= 1);
    assert (last.size() <= 1);

    unsigned entries = first.size()+last.size()+append.size()+placed.size();

    result->reserve(entries);

    result->insert(result->end(), append.begin(), append.end());

    if (first.size() == 1) result->insert(result->begin(), first[0]);
    if (last.size() == 1) result->push_back(last[0]);

    for (typename std::map<G4int, Wrapper>::iterator placedIter = placed.begin();	
	 placedIter != placed.end(); ++placedIter) {

      typename Result::iterator resultIter = result->begin() + placedIter->first;
      //      G4cout<<"jane inserting seed at "<<placedIter->first<<G4endl;
      result->insert(resultIter, placedIter->second);
    }

  }

  unsigned Size() 
  {
    return fStore.size();
  }

private:
  
  Store fStore;

};

#endif
