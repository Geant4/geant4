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
/*
 * G4KDTree.cc
 *
 *  Created on: 22 oct. 2013
 *      Author: kara
 */

#include "globals.hh"
#include <cstdio>
#include <cmath>
#include "G4KDTree.hh"
#include "G4KDMap.hh"
#include "G4KDNode.hh"
#include "G4KDTreeResult.hh"
#include <list>
#include <iostream>

using namespace std;

G4Allocator<G4KDTree>*& G4KDTree::fgAllocator()
{
  G4ThreadLocalStatic G4Allocator<G4KDTree>* _instance = nullptr;
  return _instance;
}

//______________________________________________________________________
// KDTree methods
G4KDTree::G4KDTree(size_t k)
  : fDim(k)
  ,fKDMap(new G4KDMap(k))
{}

G4KDTree::~G4KDTree()
{
  if(fRoot){
    __Clear_Rec(fRoot);
    fRoot = nullptr;
  }

  if(fRect){
    delete fRect;
    fRect = nullptr;
  }

  if(fKDMap){
    delete fKDMap;
    fKDMap = nullptr;
  }
}

void* G4KDTree::operator new(size_t)
{
  if(!fgAllocator()){
    fgAllocator() = new G4Allocator<G4KDTree>;
  }
  return (void*) fgAllocator()->MallocSingle();
}

void G4KDTree::operator delete(void* aNode)
{
  fgAllocator()->FreeSingle((G4KDTree*) aNode);
}

void G4KDTree::Print(std::ostream& out) const
{
  if(fRoot){
    fRoot->Print(out);
  }
}

void G4KDTree::Clear()
{
  __Clear_Rec(fRoot);
  fRoot    = nullptr;
  fNbNodes = 0;

  if(fRect)
  {
    delete fRect;
    fRect = nullptr;
  }
}

void G4KDTree::__Clear_Rec(G4KDNode_Base* node)
{
  if(!node)
  {
    return;
  }

  if(node->GetLeft())
  {
    __Clear_Rec(node->GetLeft());
  }
  if(node->GetRight())
  {
    __Clear_Rec(node->GetRight());
  }

  delete node;
}

void G4KDTree::__InsertMap(G4KDNode_Base* node) { fKDMap->Insert(node); }

void G4KDTree::Build()
{
  size_t Nnodes = fKDMap->GetSize();

  G4cout << "********************" << G4endl;
  G4cout << "template<typename PointT> G4KDTree<PointT>::Build" << G4endl;
  G4cout << "Map size = " << Nnodes << G4endl;

  G4KDNode_Base* root = fKDMap->PopOutMiddle(0);

  if(root == nullptr)
  {
    return;
  }

  fRoot = root;
  fNbActiveNodes++;
  fRect = new HyperRect(fDim);
  fRect->SetMinMax(*fRoot, *fRoot);

  Nnodes--;

  G4KDNode_Base* parent = fRoot;

  for(size_t n = 0; n < Nnodes; n += fDim)
  {
    for(size_t dim = 0; dim < fDim; dim++)
    {
      G4KDNode_Base* node = fKDMap->PopOutMiddle(dim);
      if(node)
      {
        parent->Insert(node);
        fNbActiveNodes++;
        fRect->Extend(*node);
        parent = node;
      }
    }
  }
}

G4KDTreeResultHandle G4KDTree::Nearest(G4KDNode_Base* node)
{
  if(!fRect)
  {
    return nullptr;
  }

  std::vector<G4KDNode_Base*> result;
  G4double dist_sq = DBL_MAX;

  /* Duplicate the bounding hyperrectangle, we will work on the copy */
  auto newrect = new HyperRect(*fRect);

  /* Search for the nearest neighbour recursively */
  G4int nbresult = 0;

  __NearestToNode(node, fRoot, *node, result, &dist_sq, newrect, nbresult);

  /* Free the copy of the hyperrect */
  delete newrect;

  /* Store the result */
  if(!result.empty())
  {
    G4KDTreeResultHandle rset(new G4KDTreeResult(this));
    G4int j = 0;
    while(j < nbresult)
    {
      rset->Insert(dist_sq, result[j]);
      j++;
    }
    rset->Rewind();

    return rset;
  }
  else
  {
    return nullptr;
  }
}

G4KDTreeResultHandle G4KDTree::NearestInRange(G4KDNode_Base* node,
                                              const G4double& range)
{
  if(!node)
  {
    return nullptr;
  }
  G4int ret(-1);

  auto* rset = new G4KDTreeResult(this);

  const G4double range_sq = sqr(range);

  if((ret = __NearestInRange(fRoot, *node, range_sq, range, *rset, 0, node)) ==
     -1)
  {
    delete rset;
    return nullptr;
  }
  rset->Sort();
  rset->Rewind();
  return rset;
}
