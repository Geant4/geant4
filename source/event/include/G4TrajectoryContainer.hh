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
// G4TrajectoryContainer
//
// Class description:
//
// This is a container of G4VTrajectory objects and the object of this
// container will be associated to G4Event object.

// Author: Makoto Asai (SLAC)
// --------------------------------------------------------------------
#ifndef G4TrajectoryContainer_hh
#define G4TrajectoryContainer_hh 1

#include <vector>

#include "globals.hh"
#include "evtdefs.hh"
#include "G4VTrajectory.hh"
#include "G4Allocator.hh"

using TrajectoryVector = std::vector<G4VTrajectory*>;

class G4TrajectoryContainer
{
  public:

    G4TrajectoryContainer();
    ~G4TrajectoryContainer();

    G4TrajectoryContainer(const G4TrajectoryContainer&) = delete;
    G4TrajectoryContainer& operator=(const G4TrajectoryContainer&) = delete;

    inline void *operator new(std::size_t);
    inline void operator delete(void* anEvent);

    G4bool operator==(const G4TrajectoryContainer& right) const;
    G4bool operator!=(const G4TrajectoryContainer& right) const;

    inline std::size_t size() const { return vect->size(); }
    inline void push_back(G4VTrajectory* p) { vect->push_back(p); }
    inline std::size_t entries() const { return size(); }
    inline G4bool insert(G4VTrajectory* p) { push_back(p); return true; }
    inline void clearAndDestroy()
    {
      for(std::size_t i=0; i<size(); ++i) delete (*vect)[i];
      vect->clear();
    }
    inline G4VTrajectory* operator[](std::size_t n) { return (*vect)[n]; }
    inline TrajectoryVector* GetVector() const { return vect; }

  private:

    TrajectoryVector* vect = nullptr;
};

extern G4EVENT_DLL
G4Allocator<G4TrajectoryContainer>*& aTrajectoryContainerAllocator();

inline void* G4TrajectoryContainer::operator new(std::size_t)
{
  if (aTrajectoryContainerAllocator() == nullptr)
  {
    aTrajectoryContainerAllocator() = new G4Allocator<G4TrajectoryContainer>;
  }
  return (void*)aTrajectoryContainerAllocator()->MallocSingle();
}

inline void G4TrajectoryContainer::operator delete(void* aTrajectoryContainer)
{
  aTrajectoryContainerAllocator()
   ->FreeSingle((G4TrajectoryContainer*)aTrajectoryContainer);
}

#endif
