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
// $Id: G4IT.hh 102616 2017-02-10 07:57:14Z gcosmo $
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

#ifndef G4IT_h
#define G4IT_h 1

#include "globals.hh"
#include "G4ITType.hh"
#include "G4ThreeVector.hh"
#include "G4VUserTrackInformation.hh"
#include "G4KDNode.hh"

///
// To implement your own IT class, you should use
// ITDef(MyClass) in the class define in your MyClass.hh
// and ITImp(MyClass) in your MyClass.cc
// For instance, see G4Molecule
///

class G4IT;
template<>
G4KDNode<G4IT>::~G4KDNode();

class G4TrackingInformation;
//template<typename PointT> class G4KDNode;
class G4KDNode_Base;
class G4ITBox;
class G4Track;

G4IT* GetIT(const G4Track* track);
G4IT* GetIT(const G4Track& track);

template<class OBJECT> class G4FastListNode;
typedef G4FastListNode<G4Track> G4TrackListNode;

/**
 * G4IT is a interface which allows the inheriting object
 * to be tracked using G4ITStepProcessor
 * The inheriting class must implement the operator < , ==
 * and != in order to enable the sorting out.
 * also the concrete header of MyIT ("MyIt.hh") should contain : ITDef(MyIT)
 * and the source of MyIT.cc : ITImp(MyIT)
 */

class G4IT : public virtual G4VUserTrackInformation
{
public:
  G4IT();
  G4IT(G4Track*);
  virtual ~G4IT();

//  inline void *operator new(size_t);
//  inline void operator delete(void *aIT);

  virtual void Print() const
  {
    ;
  }
  virtual const G4String& GetName() const = 0;

  ///
  // You should not worried of implementing diff, equal
  // and GetType.
  // When using ITDef(MyClass) this will be done.
  // However, you need to implement in the concrete class
  // even fake operators for < and ==
  // They will be used by diff and equal.
  ///
  virtual G4bool diff(const G4IT& right) const = 0;
  virtual G4bool equal(const G4IT& right) const = 0;
  G4bool operator<(const G4IT& right) const;
  G4bool operator==(const G4IT& right) const;
  G4bool operator!=(const G4IT& right) const;

  void SetTrack(G4Track*);
  inline G4Track* GetTrack();
  inline const G4Track* GetTrack() const;

  void RecordCurrentPositionNTime();
  const G4ThreeVector& GetPosition() const;
  double operator[](int i) const;
  const G4ThreeVector& GetPreStepPosition() const;
  G4double GetPreStepLocalTime() const;
  G4double GetPreStepGlobalTime() const;

  inline void SetPrevious(G4IT*);
  inline void SetNext(G4IT*);
  inline G4IT* GetPrevious();
  inline G4IT* GetNext();
  inline const G4IT* GetPrevious() const;
  inline const G4IT* GetNext() const;
  inline void SetITBox(G4ITBox*);
  inline const G4ITBox* GetITBox() const;
  void TakeOutBox();
  inline void SetNode(G4KDNode_Base*);
  inline G4KDNode_Base* GetNode() const;

  inline void SetParentID(int, int);
  inline void GetParentID(int&, int&);

  inline G4TrackingInformation* GetTrackingInfo()
  {
    return fpTrackingInformation;
  }

  inline G4TrackListNode* GetListNode()
  {
    return fpTrackNode;
  }
  inline void SetListNode(G4TrackListNode* node)
  {
    fpTrackNode = node;
  }

  virtual const G4ITType GetITType() const = 0;

  virtual G4ITType GetITSubType() const
  {
    return 0;
  }

protected:
  G4IT(const G4IT&);
  G4IT& operator=(const G4IT&);
  G4Track* fpTrack;

private:
  G4ITBox * fpITBox;
  G4IT* fpPreviousIT;
  G4IT* fpNextIT;
  G4KDNode_Base* fpKDNode;

  int fParentID_A;
  int fParentID_B;

  G4TrackingInformation* fpTrackingInformation;
  G4TrackListNode* fpTrackNode;
};
//------------------------------------------------------------------------------

inline const G4ITBox* G4IT::GetITBox() const
{
  return fpITBox;
}

inline void G4IT::SetITBox(G4ITBox * aITBox)
{
  fpITBox = aITBox;
}

inline void G4IT::SetPrevious(G4IT* aIT)
{
  fpPreviousIT = aIT;
}

inline void G4IT::SetNext(G4IT* aIT)
{
  fpNextIT = aIT;
}

inline G4IT* G4IT::GetPrevious()
{
  return fpPreviousIT;
}

inline G4IT* G4IT::GetNext()
{
  return fpNextIT;
}

inline void G4IT::SetTrack(G4Track* track)
{
  fpTrack = track;
}

inline G4Track* G4IT::GetTrack()
{
  return fpTrack;
}

inline const G4Track* G4IT::GetTrack() const
{
  return fpTrack;
}

inline void G4IT::SetParentID(int p_a, int p_b)
{
  fParentID_A = p_a;
  fParentID_B = p_b;
}

inline void G4IT::GetParentID(int& p_a, int&p_b)
{
  p_a = fParentID_A;
  p_b = fParentID_B;
}

inline const G4IT* G4IT::GetPrevious() const
{
  return fpPreviousIT;
}

inline const G4IT* G4IT::GetNext() const
{
  return fpNextIT;
}

inline void G4IT::SetNode(G4KDNode_Base* aNode)
{
  fpKDNode = aNode;
}

inline G4KDNode_Base* G4IT::GetNode() const
{
  return fpKDNode;
}
#endif
