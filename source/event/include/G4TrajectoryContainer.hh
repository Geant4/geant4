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
// $Id: G4TrajectoryContainer.hh,v 1.8 2001-07-13 15:01:47 gcosmo Exp $
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

#include "G4VTrajectory.hh"

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
      for(size_t i=0;i<size();i++) delete (*this)[i];
      clear();
    }
};
#endif
