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
// Author: Mathieu Karamitros (kara@cenbg.in2p3.fr)
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

//int RemoveNode(G4KDNode* node)
//{
//    if(!node) return 0;
//    return node->RemoveNode() ;
//}

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
    fAxis = axis;
    SetPosition(position);
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
    fTree 	= 0 ;
}

void G4KDNode::RetrieveNodeList(std::list<G4KDNode*>& output)
{
    output.push_back(this);

    if(fLeft)
        fLeft->RetrieveNodeList(output);

    if(fRight)
        fRight->RetrieveNodeList(output);
}

// NOT READY YET
//G4KDNode* G4KDNode::FindNearestOnDirection(G4KDNode* node, const int& k,
//                                           double& bestmatch, G4KDNode** bestnode)
//{
//    if(!node) return 0 ;
//    double dist = fabs(node->fPosition[k] - fPosition[k]);

//    if(dist < bestmatch)
//    {
//        bestmatch = dist ;
//        *bestnode = node ;
//    }

//    if(node->fAxis == k)
//    {
//        G4KDNode* nextNode = node->fPosition[k] < fPosition[k] ? node->fRight : node->fLeft ;
//        if(nextNode) FindNearestOnDirection(nextNode, k, bestmatch, bestnode);
//    }
//    else
//    {
//        if(node->fRight) FindNearestOnDirection(node->fRight, k, bestmatch, bestnode) ;
//        if(node->fLeft)  FindNearestOnDirection(node->fLeft,  k, bestmatch, bestnode) ;
//    }

//    return *bestnode ;
//}

//int G4KDNode::UpdateSubTree()
//{
//    G4KDNode* c_node_on_k (0) ;

//    double best_dist = DBL_MAX;

//    if(fLeft) {   FindNearestOnDirection(fLeft,  fAxis, best_dist, &c_node_on_k) ;}
//    if(fRight){   FindNearestOnDirection(fRight, fAxis, best_dist, &c_node_on_k) ;}

//    if(!fLeft && !fRight)
//    {
//        if(fSide == 1)
//        {
//            fParent->fRight = 0 ;
//        }
//        else
//            fParent->fLeft = 0 ;
//        return 0 ;
//    }
//    else if(!c_node_on_k)
//    {
//        G4String errMsg = "Something wrong : no closest nodes on direction";
//        errMsg += fAxis ;
//        errMsg += " found and yet the node has kids" ;
//        G4Exception(__PRETTY_FUNCTION__,"",FatalErrorInArgument, errMsg);
//    }
//    else
//    {
//        if(c_node_on_k->fLeft||c_node_on_k->fRight)
//        {
//            c_node_on_k->UpdateSubTree();
//        }
//        if(fSide == 1)
//        {
//            fParent->fRight = c_node_on_k ;
//        }
//        else
//            fParent->fLeft = c_node_on_k ;

//        c_node_on_k->fParent = fParent ;

//        if(fRight)
//        {
//            c_node_on_k->fRight = fRight ;
//            fRight->fParent = c_node_on_k ;
//        }
//        if(fLeft)
//        {
//            c_node_on_k->fLeft = fLeft ;
//            fLeft->fParent = c_node_on_k ;
//        }
//    }
//    return 0 ;
//}

//int G4KDNode::RemoveNode()
//{
//    if(fParent)
//    {
//        //int returned_value ;
//        UpdateSubTree();
//        fParent = 0; fLeft = 0; fRight = 0;
//        fTree = 0;
//    }
//    else
//    {
//        // Refaire un tree à partir du tree existant
//        std::list<G4KDNode*> node_lst;
//        RetrieveNodeList(node_lst);
//        PullSubTree();

//        std::list<G4KDNode*> ::iterator it = node_lst.begin();

//        G4KDNode* newRoot = *it;

//        for(it++ ; it != node_lst.end() ; it++)
//        {
//            G4KDNode* node = *it;
//            newRoot->Insert(node, node->fPosition) ;
//        }

//        fTree -> fRoot  = newRoot ;
//        this -> fLeft  = 0 ;
//        this -> fRight = 0 ;
//        this -> fTree  = 0 ;
//    }

//    return 0;
//}
