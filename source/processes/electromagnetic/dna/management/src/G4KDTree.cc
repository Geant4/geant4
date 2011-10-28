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
/*
 * This file was inspired from an existing algorithm for kdtree's
 * algorithm written by Tamas Nepusz <tamas@cs.rhul.ac.uk> and under
 * Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>
 * The original open-source version of this code
 * may be found at http://code.google.com/p/kdtree/
*/

#include "globals.hh"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <alloca.h>
#include "G4KDTree.hh"
#include "G4KDNode.hh"
#include "G4KDTreeResult.hh"
#include <list>
#include <iostream>

using namespace std;

//______________________________________________________________________
struct kdhyperrect
{
    int fDim;
    double *fMin, *fMax;              /* minimum/maximum coords */

    kdhyperrect(int dim, const double *min, const double *max)
    {
        fDim = dim;
        fMin = new double[fDim];
        fMax = new double[fDim];
        size_t size = fDim * sizeof(double);
        memcpy(fMin, min, size);
        memcpy(fMax, max, size);
   }


    ~kdhyperrect()
    {
        delete[] fMin;
        delete[] fMax;
    }

    kdhyperrect(const kdhyperrect& rect)
    {
        fDim = rect.fDim;
        fMin = new double[fDim];
        fMax = new double[fDim];
        size_t size = fDim * sizeof(double);
        memcpy(fMin, rect.fMin, size);
        memcpy(fMax, rect.fMax, size);
    }

    void Extend(const double *pos)
    {
        int i;

        for (i=0; i < fDim; i++)
        {
            if (pos[i] < fMin[i])
            {
                fMin[i] = pos[i];
            }
            if (pos[i] > fMax[i])
            {
                fMax[i] = pos[i];
            }
        }
    }

    bool CompareDistSqr(const double *pos, const double* bestmatch)
    {
        double result = 0;

        for (int i=0; i < fDim; i++)
        {
            if (pos[i] < fMin[i])
            {
                result += sqr(fMin[i] - pos[i]);
            }
            else if (pos[i] > fMax[i])
            {
                result += sqr(fMax[i] - pos[i]);
            }

            if(result >= *bestmatch) return false ;
        }

        return true ;
    }
};

//______________________________________________________________________
// KDTree methods
G4KDTree::G4KDTree (int k)
{
    fDim = k;
    fRoot = 0;
    fDestr = 0;
    fRect = 0;
    fNbNodes = 0;
}

G4KDTree::~G4KDTree ()
{
    if(fRoot) __Clear_Rec(fRoot);
    fRoot = 0;

    if (fRect)
    {
        delete fRect;
        fRect = 0;
    }
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

void G4KDTree::__Clear_Rec(G4KDNode *node)
{
    if(!node) return;

    if(node->fLeft)  __Clear_Rec(node->fLeft);
    if(node->fRight) __Clear_Rec(node->fRight);

    if(fDestr)
    {
        if(node->fData) fDestr(node->fData);
    }
    delete node;
}

G4KDNode* G4KDTree::Insert(const double *pos, void *data)
{
    G4KDNode* node = 0 ;
    if(!fRoot)
    {
        fRoot =  new G4KDNode(this,pos,data,0, 0);
        node = fRoot;
        fNbNodes = 0;
        fNbNodes++;
    }
    else
    {
        if((node=fRoot->Insert(pos, data)))
        {
            fNbNodes++;
        }
    }

    if (fRect == 0)
    {
        fRect = new kdhyperrect(fDim,pos,pos);
    }
    else
    {
        fRect->Extend(pos);
    }

    return node;
}

G4KDNode* G4KDTree::Insert(const double& x, const double& y, const double& z, void *data)
{
    double buf[3];
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
    return Insert(buf, data);
}

//__________________________________________________________________
int G4KDTree::__NearestInRange(G4KDNode *node, const double *pos, const double& range_sq,
                          const double& range, G4KDTreeResult& list, int ordered, G4KDNode *source_node)
{
    if(!node) return 0;

    double dist_sq(DBL_MAX), dx(DBL_MAX);
    int ret(-1), added_res(0);

    if(node-> fData && node != source_node)
    {
        bool do_break = false ;
        dist_sq = 0;
        for(int i=0; i<fDim; i++)
        {
            dist_sq += sqr(node->fPosition[i] - pos[i]);
            if(dist_sq > range_sq)
            {
                do_break = true;
                break;
            }
        }
        if(!do_break && dist_sq <= range_sq)
        {
            list.Insert(dist_sq, node);
            added_res = 1;
        }
    }

    dx = pos[node->fAxis] - node->fPosition[node->fAxis];

    ret = __NearestInRange(dx <= 0.0 ? node->fLeft : node->fRight, pos, range_sq, range, list, ordered, source_node);
    if(ret >= 0 && fabs(dx) <= range)
    {
        added_res += ret;
        ret = __NearestInRange(dx <= 0.0 ? node->fRight : node->fLeft, pos, range_sq, range, list, ordered, source_node);
    }

    if(ret == -1)
    {
        return -1;
    }
    added_res += ret;

    return added_res;
}

//__________________________________________________________________
void G4KDTree::__NearestToPosition(G4KDNode *node, const double *pos, G4KDNode *&result,
                         double *result_dist_sq, kdhyperrect* rect)
{
    int dir = node->fAxis;
    int i;
    double dummy(0.), dist_sq(-1.);
    G4KDNode *nearer_subtree(0), *farther_subtree (0);
    double *nearer_hyperrect_coord(0),*farther_hyperrect_coord(0);

    /* Decide whether to go left or right in the tree */
    dummy = pos[dir] - node->fPosition[dir];
    if (dummy <= 0)
    {
        nearer_subtree = node->fLeft;
        farther_subtree = node->fRight;

        nearer_hyperrect_coord = rect->fMax + dir;
        farther_hyperrect_coord = rect->fMin + dir;

        // DEBUG
//        G4cout << "rect->fMax :" << rect->fMax << G4endl;
//        G4cout << "*rect->fMax :" << *(rect->fMax) << G4endl;
//        G4cout << "rect->fMax + dir : " << rect->fMax + dir << G4endl;
//        G4cout << "*rect->fMax + dir : " << *(rect->fMax + dir) << G4endl;

//        for(int i = 0 ; i < fDim ; i++)
//        {
//            G4cout << "rect->fMax["<< i <<"] : " << rect->fMax[i] << G4endl;
//            G4cout << "rect->fMin["<< i <<"] : " << rect->fMin[i] << G4endl;
//        }

//        G4cout << "rect->fMin :" << rect->fMin << G4endl;
//        G4cout << "*rect->fMin :" << *(rect->fMin) << G4endl;
//        G4cout << "rect->fMin + dir : " << rect->fMin + dir << G4endl;
//        G4cout << "*rect->fMin + dir : " << *(rect->fMin + dir) << G4endl;
//        getchar();

    }
    else
    {
        nearer_subtree = node->fRight;
        farther_subtree = node->fLeft;
        nearer_hyperrect_coord = rect->fMin + dir;
        farther_hyperrect_coord = rect->fMax + dir;
    }

    if (nearer_subtree)
    {
        /* Slice the hyperrect to get the hyperrect of the nearer subtree */
        dummy = *nearer_hyperrect_coord;
        *nearer_hyperrect_coord = node->fPosition[dir];
        /* Recurse down into nearer subtree */
        __NearestToPosition(nearer_subtree, pos, result, result_dist_sq, rect);
        /* Undo the slice */
        *nearer_hyperrect_coord = dummy;
    }

    /* Check the distance of the point at the current node, compare it
     * with our best so far */
    if(node->fData)
    {
        dist_sq = 0;
        bool do_break = false ;
        for(i=0; i < fDim; i++)
        {
            dist_sq += sqr(node->fPosition[i] - pos[i]);
            if(dist_sq > *result_dist_sq)
            {
                do_break = true;
                break ;
            }
        }
        if (!do_break && dist_sq < *result_dist_sq)
        {
            result = node;
            *result_dist_sq = dist_sq;
        }
    }

    if (farther_subtree)
    {
        /* Get the hyperrect of the farther subtree */
        dummy = *farther_hyperrect_coord;
        *farther_hyperrect_coord = node->fPosition[dir];
        /* Check if we have to recurse down by calculating the closest
         * point of the hyperrect and see if it's closer than our
         * minimum distance in result_dist_sq. */
        if (rect->CompareDistSqr(pos,result_dist_sq))
        {
            /* Recurse down into farther subtree */
            __NearestToPosition(farther_subtree, pos, result, result_dist_sq, rect);
        }
        /* Undo the slice on the hyperrect */
        *farther_hyperrect_coord = dummy;
    }
}

G4KDTreeResultHandle G4KDTree::Nearest(const double *pos)
{
//    G4cout << "Nearest(pos)" << G4endl ;

    if (!fRect) return 0;

    G4KDNode *result(0);
    //    std::vector<G4KDNode*> result;
    double dist_sq = DBL_MAX;

    /* Duplicate the bounding hyperrectangle, we will work on the copy */
    kdhyperrect *newrect = new kdhyperrect(*fRect);

    /* Our first guesstimate is the root node */
    /* Search for the nearest neighbour recursively */
    __NearestToPosition(fRoot, pos, result, &dist_sq, newrect);

    /* Free the copy of the hyperrect */
    delete newrect;

    /* Store the result */
    if (result)
    {
        G4KDTreeResultHandle rset = new G4KDTreeResult(this);
        rset->Insert(dist_sq, result);
        rset -> Rewind();
        return rset;
    }
    else
    {
        return 0;
    }
}

//__________________________________________________________________
void G4KDTree::__NearestToNode(G4KDNode *source_node, G4KDNode *node,
                              const double *pos, std::vector<G4KDNode*>& result, double *result_dist_sq,
                              kdhyperrect* rect, int& nbresult)
{
    int dir = node->fAxis;
    double dummy, dist_sq;
    G4KDNode *nearer_subtree (0), *farther_subtree (0);
    double *nearer_hyperrect_coord (0), *farther_hyperrect_coord (0);

    /* Decide whether to go left or right in the tree */
    dummy = pos[dir] - node->fPosition[dir];
    if (dummy <= 0)
    {
        nearer_subtree = node->fLeft;
        farther_subtree = node->fRight;
        nearer_hyperrect_coord = rect->fMax + dir;
        farther_hyperrect_coord = rect->fMin + dir;
    }
    else
    {
        nearer_subtree = node->fRight;
        farther_subtree = node->fLeft;
        nearer_hyperrect_coord = rect->fMin + dir;
        farther_hyperrect_coord = rect->fMax + dir;
    }

    if (nearer_subtree)
    {
        /* Slice the hyperrect to get the hyperrect of the nearer subtree */
        dummy = *nearer_hyperrect_coord;
        *nearer_hyperrect_coord = node->fPosition[dir];
        /* Recurse down into nearer subtree */
        __NearestToNode(source_node, nearer_subtree, pos, result, result_dist_sq, rect, nbresult);
        /* Undo the slice */
        *nearer_hyperrect_coord = dummy;
    }

    /* Check the distance of the point at the current node, compare it
     * with our best so far */
    if(node-> fData && node != source_node)
    {
        dist_sq = 0;
        bool do_break = false;
        for(int i=0; i < fDim; i++)
        {
            dist_sq += sqr(node->fPosition[i] - pos[i]);
            if(dist_sq > *result_dist_sq)
            {
                do_break = true;
                break ;
            }
        }
        if(!do_break)
        {
            if (dist_sq < *result_dist_sq)
            {
                result.clear();
                nbresult = 1 ;
                result.push_back(node);
                *result_dist_sq = dist_sq;
            }
            else if(dist_sq == *result_dist_sq)
            {
                result.push_back(node);
                nbresult++;
            }
        }
    }

    if (farther_subtree)
    {
        /* Get the hyperrect of the farther subtree */
        dummy = *farther_hyperrect_coord;
        *farther_hyperrect_coord = node->fPosition[dir];
        /* Check if we have to recurse down by calculating the closest
         * point of the hyperrect and see if it's closer than our
         * minimum distance in result_dist_sq. */
        //        if (hyperrect_dist_sq(rect, pos) < *result_dist_sq)
        if (rect->CompareDistSqr(pos,result_dist_sq))
        {
            /* Recurse down into farther subtree */
            __NearestToNode(source_node, farther_subtree, pos, result, result_dist_sq, rect, nbresult);
        }
        /* Undo the slice on the hyperrect */
        *farther_hyperrect_coord = dummy;
    }
}

G4KDTreeResultHandle G4KDTree::Nearest(G4KDNode* node)
{
//    G4cout << "Nearest(node)" << G4endl ;
    if (!fRect)
    {
        G4cout << "Tree empty" << G4endl ;
        return 0;
    }

    double* pos = node->fPosition;
    std::vector<G4KDNode*> result;
    double dist_sq = DBL_MAX;

    /* Duplicate the bounding hyperrectangle, we will work on the copy */
    kdhyperrect *newrect = new kdhyperrect(*fRect);

    /* Search for the nearest neighbour recursively */
    int nbresult = 0 ;
    //    Nearest_Node_i(node, node->parent, pos, result, &dist_sq, newrect, nbresult);
    __NearestToNode(node, fRoot, pos, result, &dist_sq, newrect, nbresult);

    /* Free the copy of the hyperrect */
    delete newrect;

    /* Store the result */
    if (!result.empty())
    {
        G4KDTreeResultHandle rset(new G4KDTreeResult(this));
        int j = 0 ;
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

G4KDTreeResultHandle G4KDTree::Nearest(const double& x, const double& y, const double& z)
{
    double pos[3];
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    return Nearest(pos);
}

G4KDTreeResultHandle G4KDTree::NearestInRange(const double *pos, const double& range)
{
    int ret(-1);

    const double range_sq = sqr(range) ;

    G4KDTreeResultHandle rset = new G4KDTreeResult(this);
    if((ret = __NearestInRange(fRoot, pos, range_sq, range, *(rset()), 0)) == -1)
    {
        rset = 0;
        return rset;
    }
    rset->Sort();
    rset->Rewind();
    return rset;
}

G4KDTreeResultHandle G4KDTree::NearestInRange(const double& x,
                                                const double& y,
                                                const double& z,
                                                const double& range)
{
    double buf[3];
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
    return NearestInRange(buf, range);
}

G4KDTreeResultHandle G4KDTree::NearestInRange( G4KDNode* node, const double& range)
{
    if(!node) return 0 ;
    int ret(-1);

    G4KDTreeResult *rset = new G4KDTreeResult(this);

    const double range_sq = sqr(range) ;

    if((ret = __NearestInRange(fRoot, node->fPosition, range_sq, range, *rset, 0, node)) == -1)
    {
        delete rset;
        return 0;
    }
    rset->Sort();
    rset->Rewind();
    return rset;
}
