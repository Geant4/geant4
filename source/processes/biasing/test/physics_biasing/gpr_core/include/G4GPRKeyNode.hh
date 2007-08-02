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
// $Id: G4GPRKeyNode.hh,v 1.1 2007-08-02 18:12:06 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, August 2007. 
//
#ifndef G4GPRKEYNODE_HH
#define G4GPRKEYNODE_HH

#include "G4GPRObserverCollectionT.hh"

class G4GPRKeyNode {

public:

  G4GPRKeyNode():fState(true) {}

  void ChangeState() 
  {
    fState = !fState;
    fObserverCollection(this);
  }

  G4bool GetState() {return fState;}

  template <typename Pointer, typename PointerToMfn>
  void AddObserver(Pointer* pointer, PointerToMfn mfn) 
  {
    fObserverCollection.RegisterObserver("tmp", pointer, mfn);
  }

private:
  G4bool fState;
  
  G4GPRObserverCollectionT<G4GPRTypeList_1(G4GPRKeyNode*), G4String> fObserverCollection;
};

#endif
