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
// $Id: G4KDTreeResult.hh 101354 2016-11-15 08:27:51Z gcosmo $
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

#ifndef G4KDTREERESULT_HH
#define G4KDTREERESULT_HH

#include <list>
#include "globals.hh"
#include "G4ReferenceCountedHandle.hh"
#include "G4KDNode.hh"

class G4KDTree;
class G4KDNode_Base;
struct ResNode;
class G4KDTreeResult;

typedef G4ReferenceCountedHandle<G4KDTreeResult> G4KDTreeResultHandle;
typedef G4ReferenceCountedHandle<ResNode> ResNodeHandle;

/**
 * G4KDTreeResult enables to go through the nearest entities found
 * by G4KDTree.
 */

#define KDTR_parent std::vector<ResNode>

class G4KDTreeResult : protected KDTR_parent//protected std::list<ResNode>
{
protected:
  G4KDTree *fTree;
//  std::list<ResNode>::iterator fIterator;
  KDTR_parent::iterator fIterator;

public:
  G4KDTreeResult(G4KDTree*);
  virtual ~G4KDTreeResult();
  
  //  new/delete operators are overloded to use G4Allocator
  inline void *operator new(size_t);
#ifdef __IBMCPP__
  inline void *operator new(size_t sz, void* p)
  { return p;}
#endif
  inline void operator delete(void*);

  void Insert(double, G4KDNode_Base*);

  void Clear();

  void Sort();

  /* returns the size of the result set (in elements) */
  size_t GetSize() const;

  size_t size() const;

  /* rewinds the result set iterator */
  void Rewind();

  /* returns non-zero if the set iterator reached the end after the last element*/
  bool End();

  /* advances the result set iterator
   */
  void Next();

  /* returns the data pointer (can be null) of the current result set item
   * and optionally sets its position to the pointers(s) if not null.
   */
  template<typename PointT>
    PointT* GetItem() const;
  G4KDNode_Base* GetNode() const;
  template<typename PointT>
    PointT* GetItemNDistanceSQ(double& /*distance*/) const;
  double GetDistanceSqr() const;
};

//------------------------------------------------------------------------------
#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4ThreadLocal G4Allocator<G4KDTreeResult> *aKDTreeAllocator;
#else
extern G4DLLIMPORT G4ThreadLocal G4Allocator<G4KDTreeResult> *aKDTreeAllocator;
#endif

inline void * G4KDTreeResult::operator new(size_t)
{
  if (!aKDTreeAllocator) aKDTreeAllocator = new G4Allocator<G4KDTreeResult>;
  return (void *) aKDTreeAllocator->MallocSingle();
}

inline void G4KDTreeResult::operator delete(void * object)
{
  aKDTreeAllocator->FreeSingle((G4KDTreeResult *) object);
}
//------------------------------------------------------------------------------
template<typename PointT>
  PointT* G4KDTreeResult::GetItem() const
  {
    G4KDNode<PointT>* node = (G4KDNode<PointT>*) (GetNode());
    return node->GetPoint();
  }

template<typename PointT>
  PointT* G4KDTreeResult::GetItemNDistanceSQ(double& dist_sq) const
  {
    dist_sq = GetDistanceSqr();
    return this->GetItem<PointT>();
  }

#endif // G4KDTREERESULT_HH
