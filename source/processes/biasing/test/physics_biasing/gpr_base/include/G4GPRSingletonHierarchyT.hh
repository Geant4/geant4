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
// $Id: G4GPRSingletonHierarchyT.hh,v 1.2 2007-08-02 18:12:05 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// J. Tinslay, July 2007
//
#ifndef G4GPRSINGLETONHIERARCHYT_HH
#define G4GPRSINGLETONHIERARCHYT_HH

#include "G4GPRLinearHierarchyT.hh"

template <typename TList>
class G4GPRSingletonHierarchyT : public G4GPRLinearHierarchyT<TList> {
  
public:
  
  typedef TList TypeList;

  unsigned Size() {return G4GPRLinearHierarchyT<TList>::Size; }
  static G4GPRSingletonHierarchyT* Instance()
  {
//    G4cout<<"jane "<<fInstance<<G4endl;
    if (0 == fInstance) fInstance = new G4GPRSingletonHierarchyT();
    return fInstance;
  }

private:
  G4GPRSingletonHierarchyT(){}

  static G4GPRSingletonHierarchyT* fInstance;

};

template <class TList>
G4GPRSingletonHierarchyT<TList>* G4GPRSingletonHierarchyT<TList>::fInstance = 0;

#endif
