// -*- C++ -*-
// $Id: g4rwDetectorConstruction.hh,v 1.1 2004-12-06 11:01:14 radoone Exp $
#ifndef g4rwDetectorConstruction_H
#define g4rwDetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

#include "Saxana/SAXProcessor.h"
#include "Saxana/ProcessingConfigurator.h"

class gogdmlDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    gogdmlDetectorConstruction();
    ~gogdmlDetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();

  private:
    SAXProcessor sxp;
    ProcessingConfigurator config;
    G4VPhysicalVolume* fWorld;
};

#endif

