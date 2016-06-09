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
// $Id: G4KDNode.hh 64057 2012-10-30 15:04:49Z gcosmo $
//
// Author: Mathieu Karamitros (kara (AT) cenbg . in2p3 . fr) 
//
// WARNING : This class is released as a prototype.
// It might strongly evolve or even disapear in the next releases.
//
// History:
// -----------
// 10 Oct 2011 M.Karamitros created
//
// -------------------------------------------------------------------

#ifndef G4KDNODE_HH
#define G4KDNODE_HH

#include <list>

class G4KDTree;

/**
  * G4KDNode stores one entity in G4KDTree
  * This class is for internal use only
  */
class G4KDNode
{

public :
    // For root node :
    // parent = 0, axis = 0, side = 0
    G4KDNode(G4KDTree*, const double* /*position*/, void* /*data*/,
             G4KDNode* /*parent*/, int axis0);
    virtual ~G4KDNode();

    inline G4KDTree* GetTree();
    inline void SetTree(G4KDTree*);

    inline const double* GetPosition();

    int GetDim();

    inline int       GetAxis();
    inline void*     GetData();
    inline void      SetData(void*);
    inline G4KDNode* GetParent();
    inline G4KDNode* GetLeft();
    inline G4KDNode* GetRight();

    G4KDNode* FindParent(const double* x0);
    G4KDNode* Insert(const double* p, void* data);

    int	Insert(G4KDNode* newNode, double* p);
    int	Insert(G4KDNode* newNode, const double& x, const double& y, const double& z);
    int	Insert(G4KDNode* newNode);

    void InactiveNode();
    void PullSubTree();
    void RetrieveNodeList(std::list<G4KDNode*>& node_list);

protected :

    int SetPosition(const double* newposition);

    //°°°°°°°°°°°
    // Members
    //°°°°°°°°°°°
    double* fPosition;
    int fAxis;              // axis : x, y, z ...
    void *fData;
    int fSide ;             // left/right
    /* fSide == 0  : Is the root node
     * fSide == -1 : It is the left of the parent node
     * fSide == 1  : It is the right of the parent node
     */

    G4KDTree* fTree ;
    G4KDNode *fLeft, *fRight, *fParent;
    /* Left : fLeft->fPosition[axis] < this->fPosition[axis]
     * Right : fRight->fPosition[axis] > this->fPosition[axis]
     * Root node : fParent = 0
     */
private :
    G4KDNode(const G4KDNode& right);
    G4KDNode& operator=(const G4KDNode& right);
};

inline int G4KDNode::GetAxis()
{
    return fAxis;
}

inline void* G4KDNode::GetData()
{
    return fData;
}

inline void  G4KDNode::SetData(void* data)
{
    fData = data;
}

inline const double* G4KDNode::GetPosition()
{
    return fPosition;
}

inline G4KDNode* G4KDNode::GetParent()
{
    return fParent;
}

inline G4KDNode* G4KDNode::GetLeft()
{
    return fLeft;
}

inline G4KDNode* G4KDNode::GetRight()
{
    return fRight;
}

inline G4KDTree* G4KDNode::GetTree()
{
    return fTree;
}

inline void G4KDNode::SetTree(G4KDTree* tree)
{
    fTree = tree;
}

#endif // G4KDNODE_HH
