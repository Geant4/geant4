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
// $Id: G4TrajectoryModelMaker.hh,v 1.2 2005-10-24 14:03:36 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Template class for "factories" for making trajectory drawers.
// Class Description - End:

#ifndef G4TRAJECTORYMODELMAKER_HH
#define G4TRAJECTORYMODELMAKER_HH

#include "G4VTrajectoryModelMaker.hh"

template <class T> class G4TrajectoryModelMaker:
  public G4VTrajectoryModelMaker {

public: // With description
  
  G4TrajectoryModelMaker
  (const G4String& name,
   const G4String& commandPrefix = "/"):
    G4VTrajectoryModelMaker(name, commandPrefix)
  {}

  ~G4TrajectoryModelMaker() {}

  G4VTrajectoryModel* CreateModel();

};

template <class T>
G4VTrajectoryModel* G4TrajectoryModelMaker<T>::CreateModel()
{
  std::ostringstream oss;
  oss << fName << '-' << fID++;
  return new T(oss.str(), fCommandPrefix);
}

#endif
