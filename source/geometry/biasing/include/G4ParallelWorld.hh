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
// $Id: G4ParallelWorld.hh,v 1.6 2002-10-14 12:36:00 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4ParallelWorld
//
// Class description:
//
// Used internally by importance samplin and scoring in a "parallel"
// geometry. 
// The class relates the world volume (G4VPhysicalVolume), 
// G4VGeoDriver and G4VParallelStepper for one "parallel"
// geometry. 

// Author: Michael Dressel (Michael.Dressel@cern.ch)
// ----------------------------------------------------------------------
#ifndef G4ParallelWorld_hh
#define G4ParallelWorld_hh G4ParallelWorld_hh

#include "globals.hh"

class G4VPhysicalVolume;
class G4VParallelStepper;
class G4VPGeoDriver;

class G4ParallelWorld
{

public:  // with description

  explicit G4ParallelWorld(G4VPhysicalVolume &worldvolume);
    // initialisation and create G4ParallelStepper and a G4VPGeoDriver

  ~G4ParallelWorld();
    // delete G4ParallelStepper and a G4VPGeoDriver


  G4VParallelStepper &GetParallelStepper();
    // get the parallel stepper

  G4VPGeoDriver &GetGeoDriver();  
    // get the G4VPGeoDriver

private:

  G4ParallelWorld(const G4ParallelWorld &rhs);
    // create G4ParallelStepper and a G4VPGeoDriver
  
  G4ParallelWorld &operator=(const G4ParallelWorld &rhs);
    // create G4ParallelStepper and a G4VPGeoDriver

  G4VParallelStepper *fPstepper;
  G4VPGeoDriver *fPdriver;
};

#endif
