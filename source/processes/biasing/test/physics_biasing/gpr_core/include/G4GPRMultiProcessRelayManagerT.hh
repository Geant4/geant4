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
// $Id: G4GPRMultiProcessRelayManagerT.hh,v 1.1 2007-08-10 22:23:04 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRMULTIPROCESSRELAYMANAGERT_HH
#define G4GPRMULTIPROCESSRELAYMANAGERT_HH

#include "G4GPRPlacement.hh"
#include "G4GPRManagerT.hh"
#include "G4GPRMultiProcessRelayT.hh"
#include "G4GPRBinderFirst.hh"
#include <vector>

template <>
template <typename List>
class G4GPRManagerT< G4GPRMultiProcessRelayT<List> > {

public:

  typedef typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper SeedWrapper;
  typedef typename G4GPRProcessWrappers::Wrappers<List>::MultiProcessRelayWrapper MultiProcessRelayWrapper;
  typedef std::vector<typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper> Result;
  
  typedef std::vector< G4GPRMultiProcessRelayT<List>* > Store;


  void Register(G4GPRMultiProcessRelayT<List>* component)
  {
    fStore.push_back(component);
  }
  
  void operator()(Result*& result) 
  {
      
    typename Store::iterator iter = fStore.begin();
    G4cout<<"jane multi "<<fStore.size()<<G4endl;
    while (iter != fStore.end()) {

      if ((*iter)->IsActive()) {
	std::vector<unsigned> indices = (*iter)->GetProcessIndices();
	
	std::vector<SeedWrapper> processes;

	for (std::vector<unsigned>::iterator iterIdx = indices.begin(); iterIdx != indices.end(); ++iterIdx) {
	  processes.push_back((*result)[*iterIdx]);
	  G4cout<<"jane adding proc "<<*iterIdx<<G4endl;
	}

	G4GPRBinderFirst<MultiProcessRelayWrapper, std::vector<SeedWrapper>, typename SeedWrapper::Impl>* handle = 
	  new G4GPRBinderFirst<MultiProcessRelayWrapper, std::vector<SeedWrapper>, typename SeedWrapper::Impl>((*iter)->GetName(), (*iter)->GetWrapper(), processes);

	// Multiple process relay is to be used in place of original processes, so erase them
	for (std::vector<unsigned>::iterator iterErase = indices.begin(); iterErase != indices.end(); ++iterErase) {
	  result->erase(result->begin() + *iterErase);
	}
	
	G4cout<<"jane erased size "<<result->size()<<G4endl;
	SeedWrapper newWrapper(handle);
	G4cout<<"jane placement "<<(*iter)->Placement()<<G4endl;

	result->insert(result->begin() + (*iter)->Placement(), newWrapper);
	//		(*result)[(*iter)->Placement()] = newWrapper;
	G4cout<<"jane multi size "<<result->size()<<G4endl;
      }
      ++iter;

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
