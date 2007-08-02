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
// $Id: G4GPRKeyManagerT.hh,v 1.1 2007-08-02 18:12:06 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRKEYMANAGER_HH
#define G4GPRKEYMANAGER_HH

#include <map>
#include <vector>
#include <deque>
#include "G4GPRKeyNode.hh"

template <typename List>
class G4GPRKeyManagerT {

public:

  G4GPRKeyManagerT():fKeyChanged(true) {}

  typedef std::deque<G4bool> Key;
  typedef std::map<G4GPRKeyNode*, G4bool*> Map;

  void ChangeState(G4GPRKeyNode* node) {
    G4cout<<"jane change state"<<G4endl;
    Map::iterator iter = fMap.find(node);
    *(iter->second) = !(*(iter->second));
    
    fKeyChanged = true;
  }
  
  G4bool KeyChanged() {return fKeyChanged;}

  void ResetKeyChanged() {fKeyChanged = false;}

  void AddNode(G4GPRKeyNode* node) 
  {
    fNodeList.push_back(node);
    fKey.push_back(node->GetState());
    fMap[node] = &fKey.back();

    node->AddObserver(this, &G4GPRKeyManagerT::ChangeState);
  }
  
  const Key& GetKey() {
    return fKey;
  }

private:

  G4bool fKeyChanged;
  Key fKey;
  std::vector<G4GPRKeyNode*> fNodeList;
  Map fMap;
 
};

#endif
