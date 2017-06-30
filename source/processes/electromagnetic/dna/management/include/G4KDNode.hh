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
// $Id: G4KDNode.hh 102616 2017-02-10 07:57:14Z gcosmo $
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

#ifndef G4KDNODE_HH
#define G4KDNODE_HH

#include <list>
//#include <map>
#include <vector>
#include <deque>
#include <ostream>

// geant4
#include <G4Allocator.hh>

class G4KDTree;
class G4KDMap;

class G4KDNode_Base
{
public:
  //----------------------------
  // For root node :
  // parent = 0, axis = 0, side = 0
  G4KDNode_Base(G4KDTree*, G4KDNode_Base* /*parent*/);
  virtual ~G4KDNode_Base();
  virtual double operator[](size_t) const = 0;
  virtual void InactiveNode();
  virtual bool IsValid() const{ return true; }

  //----------------------------
  inline G4KDTree* GetTree() const {return fTree;}
  inline void SetTree(G4KDTree* tree) {fTree = tree;}

  //----------------------------
  int GetDim() const;
  inline int GetAxis() const{return fAxis;}
  inline G4KDNode_Base* GetParent(){return fParent;}
  inline G4KDNode_Base* GetLeft(){return fLeft;}
  inline G4KDNode_Base* GetRight(){return fRight;}

  //----------------------------
  template<typename Position>
    G4KDNode_Base* FindParent(const Position& x0);
    
  template<typename PointT>
    G4KDNode_Base* Insert(PointT* point);
  template<typename PointT>
    G4KDNode_Base* Insert(const PointT& point);
  int Insert(G4KDNode_Base* newNode);

  void PullSubTree();
  void RetrieveNodeList(std::list<G4KDNode_Base*>& node_list);

  void Print(std::ostream& out, int level = 0) const;

// std::vector<std::deque<G4KDNode_Base*>::iterator>*
//   GetIteratorsForSortingAlgo(G4KDMap*);
// std::map<G4KDMap*, std::vector<std::deque<G4KDNode_Base*>::iterator>>*
//   fpIteratorInSortingAlgo;

protected:
  //°°°°°°°°°°°
  // Members
  //°°°°°°°°°°°
  size_t fAxis; // axis : x, y, z ...
  int fSide; // left/right
  /* fSide == 0  : Is the root node
   * fSide == -1 : It is the left of the parent node
   * fSide == 1  : It is the right of the parent node
   */

  G4KDTree* fTree;
  G4KDNode_Base *fLeft, *fRight, *fParent;
  /* Left : fLeft->fPosition[axis] < this->fPosition[axis]
   * Right : fRight->fPosition[axis] > this->fPosition[axis]
   * Root node : fParent = 0
   */
private:
  G4KDNode_Base(const G4KDNode_Base& right);
  G4KDNode_Base& operator=(const G4KDNode_Base& right);
};

/**
 * G4KDNode stores one entity in G4KDTree
 * This class is for internal use only
 */

template<typename PointT>
  class G4KDNode : public G4KDNode_Base
  {
  public:
    //----------------------------
    // For root node :
    // parent = 0, axis = 0, side = 0
    G4KDNode(G4KDTree*, PointT* /*point*/, G4KDNode_Base* /*parent*/);
    virtual ~G4KDNode();

    void *operator new(size_t);
    void operator delete(void *);

    inline PointT* GetPoint()
    {
      return fPoint;
    }

    virtual double operator[](size_t i) const
    {
      if(fPoint == 0) abort();
      return (*fPoint)[i];
    }

    virtual void InactiveNode()
    {
      fValid = false;
      G4KDNode_Base::InactiveNode();
    }

    virtual bool IsValid() const
    {
      return fValid;
    }

  protected:
    PointT* fPoint;
    bool fValid;

  private:
    G4KDNode(const G4KDNode<PointT>& right);
    G4KDNode& operator=(const G4KDNode<PointT>& right);

    static G4ThreadLocal G4Allocator<G4KDNode<PointT>>* fgAllocator;
  };

template<typename PointT>
  G4ThreadLocal G4Allocator<G4KDNode<PointT>>*
    G4KDNode<PointT>::fgAllocator = 0;

template<typename PointT>
  void* G4KDNode<PointT>::operator new(size_t)
  {
    if(!fgAllocator) fgAllocator = new G4Allocator<G4KDNode<PointT> >;
    return (void *) fgAllocator->MallocSingle();
  }

template<typename PointT>
  void G4KDNode<PointT>::operator delete(void *aNode)
  {
    fgAllocator->FreeSingle((G4KDNode<PointT> *) aNode);
  }

/**
 * G4KDNode stores one entity in G4KDTree
 * This class is for internal use only
 */

template<typename PointCopyT>
  class G4KDNodeCopy: public G4KDNode_Base
  {
  public:
    //----------------------------
    // For root node :
    // parent = 0, axis = 0, side = 0
    G4KDNodeCopy(G4KDTree* tree,
                 const PointCopyT& point,
                 G4KDNode_Base* parent) :
     G4KDNode_Base(tree, parent)
    {
      fPoint = point;
      fValid = true;
    }

    virtual ~G4KDNodeCopy(){}

    void *operator new(size_t)
    {
      if(!fgAllocator) fgAllocator = new G4Allocator<G4KDNodeCopy<PointCopyT>>;
      return (void *) fgAllocator->MallocSingle();
    }

    void operator delete(void* aNode)
    {
      fgAllocator->FreeSingle((G4KDNodeCopy<PointCopyT>*) aNode);
    }

    inline const PointCopyT& GetPoint()
    {
      return fPoint;
    }

    virtual double operator[](size_t i) const
    {
      return fPoint[i];
    }

    virtual void InactiveNode()
    {
      fValid = false;
      G4KDNode_Base::InactiveNode();
    }

    virtual bool IsValid() const
    {
      return fValid;
    }

  protected:
    PointCopyT fPoint;
    G4bool fValid;

  private:
    G4KDNodeCopy(const G4KDNodeCopy<PointCopyT>& right) :
          G4KDNode_Base(right), fPoint(0)
    {
      fValid = false;
    }

    G4KDNodeCopy<PointCopyT>&
    operator=(const G4KDNodeCopy<PointCopyT>& right)
    {
      if(this == &right) return *this;
      fPoint = right.fPoint;
      fTree = right.fTree;
      fLeft = right.fLeft;
      fRight = right.fRight;
      fParent = right.fParent;
      fSide = right.fSide;
      fAxis = right.fAxis;
      return *this;
    }

    static G4ThreadLocal G4Allocator<G4KDNodeCopy<PointCopyT>>* fgAllocator;
  };

template<typename PointT>
  G4ThreadLocal G4Allocator<G4KDNodeCopy<PointT>>*
    G4KDNodeCopy<PointT>::fgAllocator = 0;

#include "G4KDNode.icc"

#endif // G4KDNODE_HH
