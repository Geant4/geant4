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
// $Id: G4KDTree.hh 102616 2017-02-10 07:57:14Z gcosmo $
//
// Author: Mathieu Karamitros

// The code is developed in the framework of the ESA AO7146
//
// We would be very happy hearing from you, send us your feedback! :)
//
// In order for Geant4-DNA to be maintained and still open-source,
// article citations are crucial. 
// If you use Geant4-DNA chemistry and you publish papers about your software, 
// in addition to the general paper on Geant4-DNA:
//
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// we would be very happy if you could please also cite the following
// reference papers on chemistry:
//
// J. Comput. Phys. 274 (2014) 841-882
// Prog. Nucl. Sci. Tec. 2 (2011) 503-508 

#ifndef G4KDTREE_HH
#define G4KDTREE_HH 1

#include <vector>
#include "G4KDNode.hh"
#include "G4KDTreeResult.hh"

class G4KDMap;
template<typename PointT> class G4KDNode;

//__________________________________
// Methods to act on kdnode
// Methods defined in G4KDNode.cc :
void InactiveNode(G4KDNode_Base*);
void Free(G4KDNode_Base*&);
//void* GetData(G4KDNode*);
//const double* GetNodePosition(G4KDNode_Base*);
//__________________________________

/**
 * G4KDTree is used by the ITManager to locate the neareast neighbours.
 * A kdtree sorts out node in such a way that it reduces the number of node check.
 * The results of this search can be retrieved by G4KDTreeResultHandle.
 */
class G4KDTree
{
  friend class G4KDNode_Base;
public:
  G4KDTree(size_t dim = 3);
  ~G4KDTree();
  void Clear();

  void Print(std::ostream& out = G4cout) const;
  void Build();
  void NoticeNodeDeactivation()
  {
    fNbActiveNodes--;
    if (fNbActiveNodes <= 0) Clear();
  }

  size_t GetDim() const
  {
    return fDim;
  }
  int GetNbNodes() const
  {
    return fNbNodes;
  }
  G4KDNode_Base* GetRoot()
  {
    return fRoot;
  }

  template<typename PointT>
    G4KDNode_Base* InsertMap(PointT* pos);

  // Insert and attache the data to a node at the specified position
  // In return, it gives you the corresponding node
  template<typename PointT> G4KDNode_Base* Insert(PointT* pos); // 3D

  template<typename PointT> G4KDNode_Base* Insert(const PointT& pos); // 3D

  /* Find one of the nearest nodes from the specified point.
   *
   * This function returns a pointer to a result set with at most one element.
   */
  template<typename Position> G4KDTreeResultHandle Nearest(const Position& pos);
  G4KDTreeResultHandle Nearest(G4KDNode_Base* node);

  /* Find any nearest nodes from the specified point within a range.
   *
   * This function returns a pointer to a result set, which can be manipulated
   * by the G4KDTreeResult.
   * The returned pointer can be null as an indication of an error. Otherwise
   * a valid result set is always returned which may contain 0 or more elements.
   */
  template<typename Position>
    G4KDTreeResultHandle NearestInRange(const Position& pos,
                                        const double& range);
  G4KDTreeResultHandle NearestInRange(G4KDNode_Base* node, const double& range);

  void *operator new(size_t);
  void operator delete(void *);

protected:

  //______________________________________________________________________
  class HyperRect
  {
  public:
    HyperRect(size_t dim)
    {
      fDim = dim;
      fMin = new double[fDim];
      fMax = new double[fDim];
    }

    template<typename Position>
      void SetMinMax(const Position& min, const Position& max)
      {
        for (size_t i = 0; i < fDim; i++)
        {
          fMin[i] = min[i];
          fMax[i] = max[i];
        }
      }

    ~HyperRect()
    {
      delete[] fMin;
      delete[] fMax;
    }

    HyperRect(const HyperRect& rect)
    {
      fDim = rect.fDim;
      fMin = new double[fDim];
      fMax = new double[fDim];

      for (size_t i = 0; i < fDim; i++)
      {
        fMin[i] = rect.fMin[i];
        fMax[i] = rect.fMax[i];
      }
    }

    template<typename Position>
      void Extend(const Position& pos)
      {
        for (size_t i = 0; i < fDim; i++)
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

    template<typename Position>
      bool CompareDistSqr(const Position& pos, const double* bestmatch)
      {
        double result = 0;

        for (size_t i = 0; i < fDim; i++)
        {
          if (pos[i] < fMin[i])
          {
            result += sqr(fMin[i] - pos[i]);
          }
          else if (pos[i] > fMax[i])
          {
            result += sqr(fMax[i] - pos[i]);
          }

          if (result >= *bestmatch) return false;
        }

        return true;
      }

    size_t GetDim()
    {
      return fDim;
    }
    double* GetMin()
    {
      return fMin;
    }
    double* GetMax()
    {
      return fMax;
    }

  protected:
    size_t fDim;
    double *fMin, *fMax; /* minimum/maximum coords */

  private:
    // should not be used
    HyperRect& operator=(const HyperRect& rhs)
    {
      if (this == &rhs) return *this;
      return *this;
    }
  };

protected:
  void __InsertMap(G4KDNode_Base *node);
  void __Clear_Rec(G4KDNode_Base *node);

  template<typename Position>
    int __NearestInRange(G4KDNode_Base *node,
                         const Position& pos,
                         const double& range_sq,
                         const double& range,
                         G4KDTreeResult& list,
                         int ordered,
                         G4KDNode_Base *source_node = 0);

  template<typename Position>
    void __NearestToPosition(G4KDNode_Base *node,
                             const Position& pos,
                             G4KDNode_Base *&result,
                             double *result_dist_sq,
                             HyperRect* fRect);

  template<typename Position>
    void __NearestToNode(G4KDNode_Base *source_node,
                         G4KDNode_Base *node,
                         const Position& pos,
                         std::vector<G4KDNode_Base*>& result,
                         double *result_dist_sq,
                         HyperRect* fRect,
                         int& nbresult);

protected:
  HyperRect *fRect;
  G4KDNode_Base *fRoot;
  size_t fDim;
  int fNbNodes;
  int fNbActiveNodes;
  G4KDMap* fKDMap;

  G4ThreadLocalStatic G4Allocator<G4KDTree>* fgAllocator;
};

#include "G4KDTree.icc"

#endif // G4KDTREE_HH
