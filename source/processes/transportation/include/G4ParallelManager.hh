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
// $Id: G4ParallelManager.hh,v 1.3 2002-04-10 13:14:17 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelManager
//
// Class description:
//
// Used internally by importance sampling and scoring in a "parallel"
// geometry.
// It relates G4ParallelTransport, G4ParallelWorld to one "parallel"
// geometry and one particle type.

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelManager_hh
#define G4ParallelManager_hh G4ParallelManager_hh 

#include "globals.hh"

class G4ParallelTransport;
class G4ParallelWorld;
class G4VPhysicalVolume;

class G4ParallelManager
{

public:  // with description

  G4ParallelManager(G4VPhysicalVolume &worldvolume, 
		    const G4String &particlename);
    // One G4ParallelManager for each particle type
    // and each "patallel" geometry. Create a 
    // G4ParallelWorld.

  virtual ~G4ParallelManager();
    // delete G4ParallelWorld and G4ParallelTransport if created.

  G4ParallelWorld &GetParallelWorld();
    // get the G4ParallelWorld
  G4String GetParticleName();
    // get the particle name this G4ParallelManager is responsible for
  G4ParallelTransport *CreateParallelTransport();
    // get the G4ParallelTransport
  void Initialize();
    // initialisation
private:

  G4ParallelManager(const G4ParallelManager &);
  G4ParallelManager &operator=(const G4ParallelManager &);

private:

  G4ParallelWorld *fPworld;
  G4String fParticleName;
  G4ParallelTransport *fParallelTransport;
};

#endif
