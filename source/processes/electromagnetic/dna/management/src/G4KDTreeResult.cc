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
// $Id: G4KDTreeResult.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "G4KDTreeResult.hh"
#include "G4KDNode.hh"
#include "G4KDTree.hh"

using namespace std;

struct ResNode
{
public:
    ResNode():fNode(0),fDistanceSqr(0){;}
    ResNode(double distsqr, G4KDNode* node):fNode(node),fDistanceSqr(distsqr){;}
    ResNode(const ResNode& right)
    {
        fNode = right.fNode;
        fDistanceSqr= right.fDistanceSqr;
    }
    ~ResNode(){;}

    bool operator<(const ResNode& right) const
    {
        return (fDistanceSqr < right.fDistanceSqr);
    }

    G4KDNode* GetNode() { return fNode;}
    double GetDistanceSqr() { return fDistanceSqr;}

protected:
    G4KDNode* fNode;
    double fDistanceSqr;

private:
    ResNode& operator=(const ResNode& rhs)
    {
        if(this == &rhs) return *this;
        return *this;
    }
};

// comparison
bool CompareResNode(const ResNode& left, const ResNode& right)
{
    return left < right;
}

G4KDTreeResult::G4KDTreeResult(G4KDTree* tree) : std::list<ResNode>()
{
    fTree = tree;
}

G4KDTreeResult::~G4KDTreeResult()
{
    std::list<ResNode>::erase(begin(),end());
}

void G4KDTreeResult::Insert(double pos, G4KDNode* node)
{
    std::list<ResNode>::push_back(ResNode(pos,node));
}

void G4KDTreeResult::Clear()
{
    std::list<ResNode>::erase(begin(),end());
    fIterator = std::list<ResNode>::begin();
}

void G4KDTreeResult::Sort()
{
    std::list<ResNode>::sort(CompareResNode);
}

size_t G4KDTreeResult::GetSize()
{
    return std::list<ResNode>::size();
}

size_t G4KDTreeResult::size()
{
    return std::list<ResNode>::size();
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
    fIterator++;
}

void* G4KDTreeResult::GetItem(double*& pos)
{
    if(!pos)   pos = new double[fTree->GetDim()];
    memcpy(pos, (*fIterator).GetNode()->GetPosition(), fTree->GetDim() * sizeof *pos);
    return (*fIterator).GetNode()->GetData();
}

void* G4KDTreeResult::GetItem(double& x, double& y, double& z)
{
    x = (*fIterator).GetNode()->GetPosition()[0];
    y = (*fIterator).GetNode()->GetPosition()[1];
    z = (*fIterator).GetNode()->GetPosition()[2];

    return (*fIterator).GetNode()->GetData();
}

void* G4KDTreeResult::GetItemNDistanceSQ(double& dist_sq)
{
    dist_sq = (*fIterator).GetDistanceSqr();
    return (*fIterator).GetNode()->GetData();
}

void* G4KDTreeResult::GetItemNDistanceSQ(double*& pos, double& dist_sq)
{
    dist_sq = (*fIterator).GetDistanceSqr();
    return GetItem(pos);
}

void* G4KDTreeResult::GetItemData()
{
    return (*fIterator).GetNode()->GetData();
}

double G4KDTreeResult::GetDistanceSqr()
{
    return (*fIterator).GetDistanceSqr();
}
