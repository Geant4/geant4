// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoryContainer.hh,v 1.6 2001-02-08 06:07:16 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//	G4TrajectoryContainer
//

#ifndef G4TrajectoryContainer_h
#define G4TrajectoryContainer_h 1

#include "G4VTrajectory.hh"
// G4RWTPtrOrderedVector
//#include "g4rw/tpordvec.h"

// class description:
//
//  This is a container of G4VTrajectory objects and the object of this
// container will be associated to G4Event object.

//typedef G4RWTPtrOrderedVector<G4VTrajectory> G4TrajectoryContainer;

class G4TrajectoryContainer : public G4std::vector<G4VTrajectory*>
{
  public:
    G4TrajectoryContainer() {;}
    ~G4TrajectoryContainer() {;}

  public:
    inline G4int entries() { return size(); }
    inline G4bool insert(G4VTrajectory* p) { push_back(p); return true; }
    inline void clearAndDestroy()
    {
      for(int i=0;i<size();i++)
      { delete (*this)[i]; }
      clear();
    }
};
#endif


