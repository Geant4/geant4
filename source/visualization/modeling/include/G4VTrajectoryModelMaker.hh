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
// $Id: G4VTrajectoryModelMaker.hh,v 1.1 2005-10-24 11:20:18 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Abstract base class for "factories" for making trajectory models.
// Class Description - End:

#ifndef G4VTRAJECTORYMODELMAKER_HH
#define G4VTRAJECTORYMODELMAKER_HH

#include "globals.hh"
#include "G4String.hh"

class G4VTrajectoryModel;

class G4VTrajectoryModelMaker {

public: // With description

  G4VTrajectoryModelMaker
  (const G4String& name,
   const G4String& commandPrefix):
    fName(name),
    fCommandPrefix(commandPrefix),
    fID (0)
  {}

  virtual ~G4VTrajectoryModelMaker() {}

  virtual G4VTrajectoryModel* CreateModel() = 0;

  const G4String& GetName() {return fName;}
  const G4String& GetCommandPrefix() {return fCommandPrefix;}
  G4int GetID() {return fID;}

protected:
  
  G4String fName;           // Model maker name.
  G4String fCommandPrefix;  // Model messeneger will place commands here.
  G4int fID;                // Incremented for each instance of drawer.

};

#endif
