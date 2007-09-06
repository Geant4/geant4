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
// $Id: G4GPRMaskManagerT.hh,v 1.2 2007-09-06 22:10:09 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007. 
//
#ifndef G4GPRMASKMANAGERT_HH
#define G4GPRMASKMANAGERT_HH

#include "G4GPRPlacement.hh"
#include "G4GPRManagerT.hh"
#include "G4GPRMask.hh"
#include <vector>
#include <algorithm>

template <>
class G4GPRManagerT< G4GPRMask > {

public:

  typedef std::vector< G4GPRMask* > Store;

  void Register(G4GPRMask* component)
  {
    fStore.push_back(component);
  }
  
  template <typename Result>
  void operator()(Result*& result) 
  {
    typename Store::iterator iter = fStore.begin();

    while (iter != fStore.end()) {
      
      if ((*iter)->IsActive()) {

	 std::vector<unsigned> indices = (*iter)->GetProcessIndices();
	
	// Sort into ascending order
	std::sort(indices.begin(), indices.end());

	// Erase masked out processes
	for (std::vector<unsigned>::reverse_iterator iterErase = indices.rbegin(); iterErase != indices.rend(); ++iterErase) {
	  G4cout<<"jane mask erasing "<<*iterErase<<G4endl;
	  result->erase(result->begin() + *iterErase);
	}
      }

      ++iter;
    }
      G4cout<<"jane erased size "<<result->size()<<G4endl;
  }

  unsigned Size() 
  {
    return fStore.size();
  }

private:
  
  Store fStore;

};

#endif
