//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4TrajectoryContainer.hh,v 1.9 2002-05-15 05:47:49 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4TrajectoryContainer
//
// Class description:
//
// This is a container of G4VTrajectory objects and the object of this
// container will be associated to G4Event object.
//

// ********************************************************************
#ifndef G4TrajectoryContainer_h
#define G4TrajectoryContainer_h 1

#include "globals.hh"
#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include "g4std/vector"

typedef G4std::vector<G4VTrajectory*> TrajectoryVector;

class G4TrajectoryContainer
{
  public:
    G4TrajectoryContainer();
    ~G4TrajectoryContainer();

    inline void *operator new(size_t);
    inline void operator delete(void* anEvent);

    G4int operator==(const G4TrajectoryContainer& right) const;
    G4int operator!=(const G4TrajectoryContainer& right) const;

  public:
    inline size_t size() const { return vect->size(); }
    inline void push_back(G4VTrajectory* p) { vect->push_back(p); }
    inline G4int entries() const { return size(); }
    inline G4bool insert(G4VTrajectory* p) { push_back(p); return true; }
    inline void clearAndDestroy()
    {
      for(size_t i=0;i<size();i++) delete (*vect)[i];
      vect->clear();
    }
    inline G4VTrajectory** begin() const { return vect->begin(); }
    inline G4VTrajectory** end() const { return vect->end(); }
    inline G4VTrajectory* operator[](size_t n) { return (*vect)[n]; }
    inline TrajectoryVector* GetVector() const { return vect; }

  private:
    TrajectoryVector* vect;
};


extern G4Allocator<G4TrajectoryContainer> aTrajectoryContainerAllocator;

inline void* G4TrajectoryContainer::operator new(size_t)
{
  void* aTrajectoryContainer;
  aTrajectoryContainer = (void*)aTrajectoryContainerAllocator.MallocSingle();
  return aTrajectoryContainer;
}

inline void G4TrajectoryContainer::operator delete(void* aTrajectoryContainer)
{
  aTrajectoryContainerAllocator.FreeSingle((G4TrajectoryContainer*)aTrajectoryContainer);
}


#endif
