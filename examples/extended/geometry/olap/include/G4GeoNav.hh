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
//
// $Id: G4GeoNav.hh,v 1.4 2006-06-29 17:21:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// G4GeoNav
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#ifndef G4GeoNav_h
#define G4GeoNav_h

#include "G4String.hh"
#include <vector>
#include <regex.h>

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

#include "tree.hh"

class G4VPhysicalVolume;
class OlapDetConstr;

// qt-free version of G4GeoNav

typedef Tree<G4LogicalVolume*> LVTree;

class G4GeoNav 
{

public:
   G4GeoNav( G4LogicalVolume * root );
   ~G4GeoNav();
   G4LogicalVolume * GetLV() { return theCurLV->data(); }
   
   // get next logical volume which name matches aRegExp
   //G4LogicalVolume * NextLV(const G4String & aRegExp);

   // returns a vector of LogicalVolumes which names match the aRegExp
   G4int FilterLV(const G4String & aRegExp,
                  std::vector<G4LogicalVolume*> & result, 
                  G4bool stopAtFirst=false);
   // returns a vector of LogicalVolumes reflecting the hierarchy of theCurLVItem
   G4int PathLV(std::vector<G4LogicalVolume*> & result);
   // "cd /CMS/Tr.*/.*Barrel"
   G4LogicalVolume * ChangeLV(const G4String & aRegExp);
   // get next lv as seen from the current one
   G4LogicalVolume * NextLV();
   // "pwd"
   G4int PwdLV(std::vector<G4LogicalVolume *>&);
   // "ls"
   G4int LsLV(std::vector<G4LogicalVolume *>&);
   // theSubreeLVItem = theCurLVItem
   //void CurLVToSubtree() {theSubtreeLVItem = theCurLVItem;};
   //OlapDetConstr * GetDetConstr();
   const LVTree::node_t & root() const { return *(theLVTree->root()); }
   
protected:
   
   void RecursiveFill(LVTree::node_t*);
   
   void CountSubtreeNodes(LVTree::node_t*,G4int&);
     // recursive counting of nodes

   G4int Tokenize(const G4String &, std::vector<G4String>&);

   void FindLV(regex_t *, LVTree::node_t * ,
               std::vector<G4LogicalVolume*>&, G4bool stopAtFirst=false) ;

   void populateLVTree(LVTree::node_t*);  // recursive population of LV-Tree

   LVTree * theLVTree ;
   LVTree::node_t * theCurLV;  // used in the NewWorld
   LVTree::node_t * theRootLV; // root of the 'whole' geometry
   LVTree::node_t * theSubtreeEndLVItem; // (un)used ...
};

#endif 
