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
// $Id: G4ParallelWorld.hh,v 1.2 2002-04-09 16:23:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelWorld
//
// Class description:
//
// <<insert the description here>>

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelWorld_hh
#define G4ParallelWorld_hh

#include "globals.hh"

class G4VPhysicalVolume;
class G4VParallelStepper;
class G4VPGeoDriver;

class G4ParallelWorld
{

public:  // with description

  G4ParallelWorld(G4VPhysicalVolume &worldvolume);
  ~G4ParallelWorld();

  G4ParallelWorld(const G4ParallelWorld &rhs);
  G4ParallelWorld &operator=(const G4ParallelWorld &rhs);

  G4VPhysicalVolume *GetWorldVolume() const;
  G4VParallelStepper &GetParallelStepper();
  G4VPGeoDriver &GetGeoDriver();  

private:

  G4VPhysicalVolume *fWorldVolume;
  G4VParallelStepper *fPstepper;
  G4VPGeoDriver *fPdriver;
};

#endif
