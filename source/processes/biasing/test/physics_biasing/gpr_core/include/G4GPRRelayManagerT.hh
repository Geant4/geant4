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
// $Id: G4GPRRelayManagerT.hh,v 1.1 2007-08-08 20:50:55 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRRELAYMANAGERT_HH
#define G4GPRRELAYMANAGERT_HH

#include "G4GPRPlacement.hh"
#include "G4GPRManagerT.hh"
#include "G4GPRRelayT.hh"
#include "G4GPRBinderFirst.hh"
#include <vector>

template <>
template <typename List>
class G4GPRManagerT< G4GPRRelayT<List> > {

public:

  typedef typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper SeedWrapper;
  typedef typename G4GPRProcessWrappers::Wrappers<List>::RelayWrapper RelayWrapper;
  typedef std::vector<typename G4GPRProcessWrappers::Wrappers<List>::SeedWrapper> Result;
  
  typedef std::vector< G4GPRRelayT<List>* > Store;


  void Register(G4GPRRelayT<List>* component)
  {
    fStore.push_back(component);
  }
  
  void operator()(Result*& result) 
  {
    typename Store::iterator iter = fStore.begin();

    while (iter != fStore.end()) {

      if ((*iter)->IsActive()) {
	G4int idx = (*iter)->Placement();
	
	G4GPRBinderFirst<RelayWrapper>* handle = new G4GPRBinderFirst<RelayWrapper>((*iter)->GetName(), (*iter)->GetWrapper(), (*result)[idx]);
	
	SeedWrapper newWrapper(handle);

	(*result)[idx] = newWrapper;
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
