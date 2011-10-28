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

//*********************************************
// This class is only for internal use
class G4KDNode
{
    friend class G4KDTree ;

protected :
    double* fPosition;
    int fAxis; // axis : x, y, z ...
    void *fData;
    int fSide ; // left/right
    G4KDTree* fTree ;
    G4KDNode *fLeft, *fRight, *fParent;	/* negative/positive side */

    int SetPosition(const double* newposition);

public :
    // For root node :
    // parent = 0, axis = 0
    G4KDNode(G4KDTree*, const double* /*position*/, void* /*data*/,
             G4KDNode* /*parent*/, int axis0);
    virtual ~G4KDNode();

    G4KDTree* GetTree();
    void SetTree(G4KDTree*);

    inline const double* GetPosition();

    int GetDim();

    inline void* GetData();
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

    //    G4KDNode* FindNearestOnDirection(G4KDNode* node, const int& k, double& bestmatch, G4KDNode** bestnode);
    //    int  RemoveNode();
    //    int  UpdateSubTree();
};

inline void* G4KDNode::GetData()
{
    return fData;
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

#endif // G4KDNODE_HH
