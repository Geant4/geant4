// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TRTDetectorConstruction.hh,v 1.1 2000-07-14 15:51:18 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef TRTDetectorConstruction_h
#define TRTDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "Randomize.hh"

class TRTTrackerSD;
class SCTTrackerSD;
class HepJamesRandom;

class TRTDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    TRTDetectorConstruction();
    ~TRTDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

  private:
    TRTTrackerSD *TRTStrawSD;
    SCTTrackerSD *SCTStripSD;
    HepJamesRandom theJamesEngine;    
};

#endif
