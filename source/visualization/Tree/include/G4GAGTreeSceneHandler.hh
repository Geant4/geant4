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
// Satoshi Tanaka  31th April 2001
// A scene handler to dump geometry hierarchy GAG


#ifndef G4GAGTREESCENEHANDLER_HH
#define G4GAGTREESCENEHANDLER_HH

#include "G4VTreeSceneHandler.hh"

#include <set>

class G4VPhysicalVolume;

/////////////////////////////
#include "G4GAGTreeList.hh"
class G4String ;
/////////////////////////////

///////////////////////////////////////////////////////////////////
// The algorithm of re-generating a geometry tree 
//  of physical volumes (PV):
// 
// In each addition of the N-th PV with a given depth D, P(N;D) to
// the tree, we use the following rules:
// 
//   R1. We can add P(N;D) to one of the previously-generated PV nodes,
//       P(N_mother <N; D_mother), if D_mother = D - 1.  
// 
//   R2. The mother PV node P(N_mother, D_mother) of P(N,D)
//       is searched in the direction to the root node 
//       (the world volume), starting from P(N-1;*).
// 
//   R3. The first PV node found in R2 is adopted as the mother.
// 
// This rule can be easily implemented with a stack.
// 
///////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// The output to G4cout becomes, for example, 
//
//  @@DTREEBegin
//  /WORLD.0/PV1.1 
//  /WORLD.0/PV1.1/PV2.2
//  /WORLD.0/PV1.1/PV3.3
//  /WORLD.0/PV1.1/PV3.3/PV4.4
//  /WORLD.0/PV5.5
//  ...
//  @@DTREEBegin
//
//  where ".number" is the index label added by the 
//  G4GAGTree driver.
///////////////////////////////////////////////////////////////////

class G4GAGTreeSceneHandler: public G4VTreeSceneHandler {
public:
  G4GAGTreeSceneHandler(G4VGraphicsSystem& system,
			  const G4String& name);
  virtual ~G4GAGTreeSceneHandler();
  void BeginModeling();
  void EndModeling();
protected:
  void RequestPrimitives(const G4VSolid&);
  // Overrides G4VScenehandler::RequestPrimitives and implements dump
  // of the geometry hierarchy.
  std::set<G4LogicalVolume*,std::less<G4LogicalVolume*> > fLVSet;
  typedef
  std::set<G4LogicalVolume*,std::less<G4LogicalVolume*> >::iterator
  LVSetIterator;
  std::set<G4VPhysicalVolume*,std::less<G4VPhysicalVolume*> > fReplicaSet;
  typedef
  std::set<G4VPhysicalVolume*,std::less<G4VPhysicalVolume*> >::iterator
  ReplicaSetIterator;

/////////////////////////////
  void ClearPVList(void);
  void InitializePVList(void);
  G4GAGTreeList<G4String>  fPVNameList ;
  G4int                    fPrevDepth  ;
  G4String                 fPrevAbsPVName ;
  G4int                    fPVCounter ;
/////////////////////////////
};

#endif
