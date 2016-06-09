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
//
// $Id: G4ParallelWorld.hh,v 1.7 2006/06/29 18:16:27 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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
