// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: gogdmlDetectorConstruction.hh,v 1.1.1.1 2002-05-31 00:34:43 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef gogdmlDetectorConstruction_H
#define gogdmlDetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

#include "SAXProcessor.hh"
#include "ProcessingConfigurator.hh"

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

