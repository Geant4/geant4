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
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "G4KDTree.hh"
#include "G4KDMap.hh"
#include "G4KDNode.hh"
#include "G4KDTreeResult.hh"
#include <list>
#include <iostream>

using namespace std;

G4ThreadLocal G4Allocator<G4KDTree>* G4KDTree::fgAllocator(0);

//______________________________________________________________________
// KDTree methods
G4KDTree::G4KDTree(size_t k) :
    fKDMap(new G4KDMap(k))
{
  fDim = k;
  fRoot = 0;
  fRect = 0;
  fNbNodes = 0;
  fNbActiveNodes = 0;
}

G4KDTree::~G4KDTree()
{
  if (fRoot) __Clear_Rec(fRoot);
  fRoot = 0;

  if (fRect)
  {
    delete fRect;
    fRect = 0;
  }

  if (fKDMap) delete fKDMap;
}

void* G4KDTree::operator new(size_t)
{
  if (!fgAllocator) fgAllocator = new G4Allocator<G4KDTree>;
  return (void *) fgAllocator->MallocSingle();
}

void G4KDTree::operator delete(void *aNode)
{
  fgAllocator->FreeSingle((G4KDTree*) aNode);
}

void G4KDTree::Print(std::ostream& out) const
{
  if (fRoot) fRoot->Print(out);
}

void G4KDTree::Clear()
{
  __Clear_Rec(fRoot);
  fRoot = 0;
  fNbNodes = 0;

  if (fRect)
  {
    delete fRect;
    fRect = 0;
  }
}

void G4KDTree::__Clear_Rec(G4KDNode_Base* node)
{
  if (!node) return;

  if (node->GetLeft()) __Clear_Rec(node->GetLeft());
  if (node->GetRight()) __Clear_Rec(node->GetRight());

  delete node;
}

void G4KDTree::__InsertMap(G4KDNode_Base *node)
{
  fKDMap->Insert(node);
}

void G4KDTree::Build()
{
  size_t Nnodes = fKDMap->GetSize();

  G4cout << "********************" << G4endl;
  G4cout << "template<typename PointT> G4KDTree<PointT>::Build" << G4endl;
  G4cout << "Map size = " << Nnodes << G4endl;

  G4KDNode_Base* root = fKDMap->PopOutMiddle(0);

  if(root == 0) return;

  fRoot = root;
  fNbActiveNodes++;
  fRect = new HyperRect(fDim);
  fRect->SetMinMax(*fRoot, *fRoot);

  Nnodes--;

  G4KDNode_Base* parent = fRoot;

  for (size_t n = 0; n < Nnodes; n += fDim)
  {
    for (size_t dim = 0; dim < fDim; dim++)
    {
      G4KDNode_Base* node = fKDMap->PopOutMiddle(dim);
      if (node)
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
  //    G4cout << "Nearest(node)" << G4endl ;
  if (!fRect)
  {
    G4cout << "Tree empty" << G4endl;
    return 0;
  }

  std::vector<G4KDNode_Base* > result;
  double dist_sq = DBL_MAX;

  /* Duplicate the bounding hyperrectangle, we will work on the copy */
  HyperRect *newrect = new HyperRect(*fRect);

  /* Search for the nearest neighbour recursively */
  int nbresult = 0;

  __NearestToNode(node, fRoot, *node, result, &dist_sq, newrect, nbresult);

  /* Free the copy of the hyperrect */
  delete newrect;

  /* Store the result */
  if (!result.empty())
  {
    G4KDTreeResultHandle rset(new G4KDTreeResult(this));
    int j = 0;
    while (j<nbresult)
    {
      rset->Insert(dist_sq, result[j]);
      j++;
    }
    rset->Rewind();

    return rset;
  }
  else
  {
    return 0;
  }
}

G4KDTreeResultHandle G4KDTree::NearestInRange(G4KDNode_Base* node,
                                              const double& range)
{
  if (!node) return 0;
  int ret(-1);

  G4KDTreeResult *rset = new G4KDTreeResult(this);

  const double range_sq = sqr(range);

  if ((ret = __NearestInRange(fRoot, *node, range_sq, range, *rset, 0, node)) == -1)
  {
    delete rset;
    return 0;
  }
  rset->Sort();
  rset->Rewind();
  return rset;
}
