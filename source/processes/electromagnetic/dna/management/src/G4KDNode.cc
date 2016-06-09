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
// $Id: G4KDNode.cc 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#include "globals.hh"
#include "G4KDNode.hh"
#include "G4KDTree.hh"

//*********************************************

//______________________________________________________________________
// Node functions

void* GetData(G4KDNode* node)
{
    return node->GetData() ;
}

const double* GetNodePosition(G4KDNode* node)
{
    return node->GetPosition() ;
}

//______________________________________________________________________

void InactiveNode(G4KDNode* node)
{
    if(!node) return ;
    node->InactiveNode() ;
}

void Free(G4KDNode*& node)
{
    if(node)
    delete node ;
    node = 0;
}
//*********************************************

G4KDNode::G4KDNode(G4KDTree* tree, const double* position, void* data,
                   G4KDNode* parent, int axis) :
    fPosition(0), fData(data), fTree(tree),
    fLeft(0), fRight(0), fParent(parent)
{
    fSide = 0;
    fAxis = axis;
    SetPosition(position);
}

// Copy constructor should not be used
G4KDNode::G4KDNode(const G4KDNode& ):
    fPosition(0), fData(0), fTree(0),
    fLeft(0), fRight(0), fParent(0)
{
    fSide = 0;
    fAxis = 0;
    fPosition = 0;
}

// Assignement should not be used
G4KDNode& G4KDNode::operator=(const G4KDNode& right)
{
    if (this == &right) return *this;
    fPosition = 0;
    fData = right.fData;
    fTree = right.fTree;
    fLeft = right.fLeft;
    fRight = right.fRight;
    fParent = right.fParent;
    fSide = right.fSide;
    fAxis = right.fAxis;
    SetPosition(right.fPosition);
    return *this;
}

G4KDNode::~G4KDNode()
{
    delete[] fPosition;
}

int G4KDNode::GetDim()
{
    if(fTree)
        return fTree->GetDim();
    else
        return -1;
}

void G4KDNode::InactiveNode()
{
    fData = 0 ;
}

int G4KDNode::SetPosition(const double* newposition)
{
    if(!newposition) return -1;
    if(!fPosition)
    {
        fPosition = new double[fTree->fDim];
    }

    memcpy(fPosition, newposition, fTree->fDim * sizeof(double));

    return 0;
}

G4KDNode* G4KDNode::FindParent(const double* x0)
{
    G4KDNode* aParent = 0 ;
    G4KDNode* next = this ;
    int split = -1 ;
    while(next)
    {
        split = next->fAxis  ;
        aParent = next ;

        if(x0[split] > next->fPosition[split])
            next = next->fRight ;
        else
            next = next->fLeft ;
    }
    return aParent ;
}

G4KDNode* G4KDNode::Insert(const double* p, void* data)
{
    G4KDNode* aParent = FindParent(p);
    // TODO check p == aParent->pos
    // Exception

    G4KDNode* newNode = new G4KDNode(fTree, p, data, aParent,
                                     aParent->fAxis +1 < fTree->fDim? aParent->fAxis+1:0);

    if(p[aParent->fAxis] > aParent->fPosition[aParent->fAxis])
    {
        aParent->fRight = newNode ;
        newNode->fSide = 1 ;
    }
    else
    {
        aParent->fLeft = newNode ;
        newNode->fSide = -1 ;
    }

    return newNode ;
}


int G4KDNode::Insert(G4KDNode* newNode, double* p)
{
    G4KDNode* aParent = FindParent(p);
    // TODO check p == aParent->pos
    // Exception

    newNode->fAxis = aParent->fAxis +1 < fTree->fDim? aParent->fAxis+1:0;
    newNode->fParent = aParent ;

    if(p[aParent->fAxis] > aParent->fPosition[aParent->fAxis])
    {
        aParent->fRight = newNode ;
        newNode->fSide = 1 ;
    }
    else
    {
        aParent->fLeft = newNode ;
        newNode->fSide = -1 ;
    }

    newNode->fRight = 0;
    newNode->fLeft = 0;

    return 0 ;
}

int G4KDNode::Insert(G4KDNode* newNode, const double& x, const double& y, const double& z)
{
    double p[3] ;
    p[0] = x;
    p[1] = y ;
    p[2] = z ;
    return Insert(newNode, p);
}

int G4KDNode::Insert(G4KDNode* newNode)
{
    return Insert(newNode, newNode->fPosition);
}

void G4KDNode::PullSubTree()
{
    if(fParent)
    {
        if(fSide == -1)
        {
            fParent->fLeft = 0;
        }
        else
            fParent->fRight = 0;
    }
    if(fLeft) fLeft -> PullSubTree();
    if(fRight) fRight-> PullSubTree();

    fParent  = 0 ;
    fRight   = 0 ;
    fLeft    = 0 ;
    fTree    = 0 ;
}

void G4KDNode::RetrieveNodeList(std::list<G4KDNode*>& output)
{
    output.push_back(this);

    if(fLeft)
        fLeft->RetrieveNodeList(output);

    if(fRight)
        fRight->RetrieveNodeList(output);
}
