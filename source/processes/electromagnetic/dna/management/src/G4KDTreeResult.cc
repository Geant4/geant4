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
// $Id: G4KDTreeResult.cc 101354 2016-11-15 08:27:51Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4KDTreeResult.hh"

using namespace std;

G4ThreadLocal G4Allocator<G4KDTreeResult> *aKDTreeAllocator= 0;

struct ResNode
{
public:
  ResNode():fNode(0),fDistanceSqr(0){;}
  ResNode(double distsqr, G4KDNode_Base* node):
  fNode(node),fDistanceSqr(distsqr)
  {;}
  
  ResNode(const ResNode& right)
  {
    fNode = right.fNode;
    fDistanceSqr= right.fDistanceSqr;
  }
  ResNode& operator=(const ResNode& rhs)
  {
    if(this == &rhs) return *this;
    fNode = rhs.fNode;
    fDistanceSqr= rhs.fDistanceSqr;
    return *this;
  }
  ~ResNode(){;}
  
  bool operator<(const ResNode& right) const
  {
    return (fDistanceSqr < right.fDistanceSqr);
  }
  
  G4KDNode_Base* GetNode() { return fNode;}
  double GetDistanceSqr() { return fDistanceSqr;}
  
protected:
  G4KDNode_Base* fNode;
  double fDistanceSqr;
};

// comparison
bool CompareResNode(const ResNode& left,
                    const ResNode& right)
{
    return left < right;
}

G4KDTreeResult::G4KDTreeResult(G4KDTree* tree) :
//std::list<ResNode>()
KDTR_parent()
{
    fTree = tree;
}

G4KDTreeResult::~G4KDTreeResult()
{
  KDTR_parent::erase(begin(),end());
    //std::list<ResNode>::erase(begin(),end());
}

void G4KDTreeResult::Insert(double dis_sq, G4KDNode_Base* node)
{
    //std::list<ResNode>::push_back(ResNode(dis_sq,node));
 KDTR_parent::push_back(ResNode(dis_sq,node));
}

void G4KDTreeResult::Clear()
{
//    std::list<ResNode>::erase(begin(),end());
//    fIterator = std::list<ResNode>::begin();
  KDTR_parent::erase(begin(),end());
  fIterator = KDTR_parent::begin();
}

void G4KDTreeResult::Sort()
{
    //std::list<ResNode>::sort(CompareResNode);
  std::sort(begin(), end(), CompareResNode);
}

size_t G4KDTreeResult::GetSize() const
{
    //return std::list<ResNode>::size();
  return KDTR_parent::size();
}

size_t G4KDTreeResult::size() const
{
  return KDTR_parent::size();
  //return std::list<ResNode>::size();
}

void G4KDTreeResult::Rewind()
{
  fIterator = begin();
}

bool G4KDTreeResult::End()
{
    return (fIterator == end());
}

void G4KDTreeResult::Next()
{
    ++fIterator;
}

double G4KDTreeResult::GetDistanceSqr() const
{
    return (*fIterator).GetDistanceSqr();
}

G4KDNode_Base* G4KDTreeResult::GetNode() const {
	return (*fIterator).GetNode();
}
