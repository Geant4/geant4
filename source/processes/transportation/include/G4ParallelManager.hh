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
// $Id: G4ParallelManager.hh,v 1.2 2002-04-09 17:40:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelManager
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelManager_hh
#define G4ParallelManager_hh

#include "globals.hh"

class G4ParallelTransport;
class G4ParallelWorld;
class G4VPhysicalVolume;

class G4ParallelManager
{

public:  // with description

  G4ParallelManager(G4VPhysicalVolume &worldvolume, 
		    const G4String &particlename);
  virtual ~G4ParallelManager();

  G4ParallelWorld &GetParallelWorld();
  G4String GetParticleName();
  G4ParallelTransport *CreateParallelTransport();
  void Initialize();

private:

  G4ParallelManager(const G4ParallelManager &);
  G4ParallelManager &operator=(const G4ParallelManager &);

private:

  G4ParallelWorld *fPworld;
  G4String fParticleName;
  G4ParallelTransport *fParallelTransport;
};

#endif
