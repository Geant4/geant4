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
// $Id: G4IT.hh 65022 2012-11-12 16:43:12Z gcosmo $
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

#ifndef G4IT_h
#define G4IT_h 1

#include "globals.hh"
#include "G4ITType.hh"
#include "G4ThreeVector.hh"
#include "G4VUserTrackInformation.hh"
#include "G4TrackingInformation.hh"

///
// To implement your own IT class, you should use
// ITDef(MyClass) in the class define in your MyClass.hh
// and ITImp(MyClass) in your MyClass.cc
// For instance, see G4Molecule
///

class G4IT;
class G4KDNode;
class G4ITBox;
class G4Track;

G4IT* GetIT(const G4Track* track) ;
G4IT* GetIT(const G4Track& track) ;

#if defined G4EM_ALLOC_EXPORT
extern G4DLLEXPORT G4Allocator<G4IT> aITAllocator;
#else
extern G4DLLIMPORT G4Allocator<G4IT> aITAllocator;
#endif

class G4TrackListNode;

/**
  * G4IT is a interface which allows the inheriting object :
  * - to be included in ITManager for the search of nearest
  * neighbour
  * - to be tracked using G4ITStepManager
  * The inheriting class must implement the operator < , ==
  * and != in order to enable the sorting out.
  * also the concrete header of MyIT ("MyIt.hh") should contain : ITDef(MyIT)
  * and the source of MyIT.cc : ITImp(MyIT)
  */

class G4IT : public virtual G4VUserTrackInformation
{
public :
    G4IT();
    G4IT(G4Track*);
    virtual ~G4IT();

    inline void *operator new(size_t);
    inline void operator delete(void *aIT);

    virtual void Print() const {}
    virtual const G4String& GetName() const = 0 ;

    ///
    // You should not worried of implementing diff, equal
    // and GetType.
    // When using ITDef(MyClass) this will be done.
    // However, you need to implement in the concrete class
    // even fake operators for < and ==
    // They will be used by diff and equal.
    ///
    virtual G4bool diff(const G4IT& right) const  = 0 ;
    virtual G4bool equal(const G4IT& right) const  = 0 ;
    G4bool operator<(const G4IT& right) const;
    G4bool operator==(const G4IT& right) const;
    G4bool operator!=(const G4IT& right) const;

    void SetTrack(G4Track*);
    inline G4Track* GetTrack();
    inline const G4Track* GetTrack() const;

    void RecordCurrentPositionNTime();

    inline void SetPrevious(G4IT*);
    inline void SetNext(G4IT*);
    inline G4IT* GetPrevious();
    inline G4IT* GetNext();
    inline const G4IT* GetPrevious() const;
    inline const G4IT* GetNext() const;
    inline void SetITBox(G4ITBox*);
    inline const G4ITBox* GetITBox() const;
    void TakeOutBox();
    inline void SetNode(G4KDNode*);

    inline void SetParentID(int,int);
    inline void GetParentID(int&,int&);

    inline const G4ThreeVector&    GetPreStepPosition() const;
    inline G4double                GetPreStepLocalTime() const;
    inline G4double                GetPreStepGlobalTime() const;
    inline G4KDNode*               GetNode() const;

    inline G4TrackingInformation* GetTrackingInfo(){return &fTrackingInformation;}

    inline G4TrackListNode* GetTrackListNode(){return fpTrackNode;}
    inline void SetTrackListNode(G4TrackListNode* node){ fpTrackNode = node;}

    virtual const G4ITType GetITType() const = 0 ;

protected :
    G4IT(const G4IT&);
    G4IT& operator=(const G4IT&);
    G4Track* fpTrack ;

private :
    G4ITBox *   fpITBox;
    G4IT*       fpPreviousIT;
    G4IT*       fpNextIT;
    G4KDNode*   fpKDNode ;

    int fParentID_A;
    int fParentID_B;

    G4TrackingInformation fTrackingInformation ;
    G4TrackListNode* fpTrackNode;
};
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
///
// Inline methods
///
inline void* G4IT::operator new(size_t)
{
    void *aIT;
    aIT = (void *) aITAllocator.MallocSingle();
    return aIT;
}

inline void G4IT::operator delete(void *aIT)
{ aITAllocator.FreeSingle((G4IT *) aIT);}

inline const G4ITBox* G4IT::GetITBox() const
{
    return fpITBox ;
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

inline void G4IT::GetParentID(int& p_a,int&p_b)
{
    p_a = fParentID_A;
    p_b = fParentID_B ;
}

inline G4double G4IT::GetPreStepGlobalTime() const
{
    return fTrackingInformation.GetPreStepGlobalTime();
}

inline G4double G4IT::GetPreStepLocalTime() const
{
    return fTrackingInformation.GetPreStepLocalTime();
}

inline const G4ThreeVector& G4IT::GetPreStepPosition() const
{
    return fTrackingInformation.GetPreStepPosition();
}

inline const G4IT* G4IT::GetPrevious() const
{
    return fpPreviousIT ;
}

inline const G4IT* G4IT::GetNext() const
{
    return fpNextIT ;
}

inline void G4IT::SetNode(G4KDNode* aNode)
{
    fpKDNode = aNode ;
}

inline G4KDNode* G4IT::GetNode() const
{
    return fpKDNode ;
}
#endif



