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
// $Id: G4KDTree.hh 64057 2012-10-30 15:04:49Z gcosmo $
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

#ifndef G4KDTREE_HH
#define G4KDTREE_HH

#include <vector>
#include "G4KDTreeResult.hh"

//__________________________________
// Methods to act on kdnode
// Methods defined in G4KDNode.cc :
void InactiveNode(G4KDNode*);
void Free(G4KDNode*&);
void* GetData(G4KDNode*);
const double* GetNodePosition(G4KDNode*);
//__________________________________

/**
  * G4KDTree is used by the ITManager to locate the neareast neighbours.
  * A kdtree sorts out node in such a way that it reduces the number of node check.
  * The results of this search can be retrieved by G4KDTreeResultHandle.
  */

class G4KDTree
{
    friend class G4KDNode ;
    int fDim;
    struct HyperRect *fRect;
    void (*fDestr)(void*);
    int fNbNodes;

protected :
    G4KDNode *fRoot;

public :
    G4KDTree(int dim = 3);
    virtual ~G4KDTree();

    void Clear();

    inline int GetDim();
    inline void SetDataDestructor(void (*fDestr)(void*));

    int GetNbNodes()    { return fNbNodes;  }
    G4KDNode* GetRoot() { return fRoot ;    }

    // Insert and attache the data to a node at the specified position
    // In return, it gives you the corresponding node
    G4KDNode* Insert(const double *pos, void *data);
    G4KDNode* Insert(const double& x, const double& y, const double& z, void *data); // 3D

    /* Find one of the nearest nodes from the specified point.
     *
     * This function returns a pointer to a result set with at most one element.
     */
    G4KDTreeResultHandle Nearest( const double *pos);
    G4KDTreeResultHandle Nearest( const double& x, const double& y, const double& z); // 3D
    G4KDTreeResultHandle Nearest( G4KDNode* node);

    /* Find any nearest nodes from the specified point within a range.
     *
     * This function returns a pointer to a result set, which can be manipulated
     * by the G4KDTreeResult.
     * The returned pointer can be null as an indication of an error. Otherwise
     * a valid result set is always returned which may contain 0 or more elements.
     */
    G4KDTreeResultHandle NearestInRange( const double *pos, const double& range);
    G4KDTreeResultHandle NearestInRange( const double& x,
                                    const double& y,
                                    const double& z,
                                    const double& range); // 3D
    G4KDTreeResultHandle NearestInRange( G4KDNode* node, const double& range);

protected :
    void __Clear_Rec(G4KDNode *node) ;

    int __NearestInRange(G4KDNode *node,
                         const double *pos,
                         const double& range_sq,
                         const double& range,
                         G4KDTreeResult& list,
                         int ordered,
                         G4KDNode *source_node = 0);

    void __NearestToPosition(G4KDNode *node,
                             const double *pos,
                             G4KDNode *&result,
                             double *result_dist_sq,
                             struct HyperRect* fRect);

    void __NearestToNode(G4KDNode *source_node,
                         G4KDNode *node,
                         const double *pos,
                         std::vector<G4KDNode*>& result,
                         double *result_dist_sq,
                         struct HyperRect* fRect,
                         int& nbresult) ;
};

inline int G4KDTree::GetDim()
{
    return fDim ;
}

void G4KDTree::SetDataDestructor(void (*fct)(void*))
{
    fDestr = fct;
}

#endif // G4KDTREE_HH
