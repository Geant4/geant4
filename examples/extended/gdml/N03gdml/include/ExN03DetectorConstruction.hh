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
// $Id: ExN03DetectorConstruction.hh,v 1.2 2002-06-03 12:09:31 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// Comments
//
// --------------------------------------------------------------

#ifndef ExN03DetectorConstruction_H
#define ExN03DetectorConstruction_H 1

#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"

#include "SAXProcessor.hh"
#include "ProcessingConfigurator.hh"

class ExN03CalorimeterSD;
class ExN03DetectorMessenger;

class ExN03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN03DetectorConstruction();
    ~ExN03DetectorConstruction();

  public:
         G4VPhysicalVolume* Construct();
  
         G4LogicalVolume*   FindLogicalVolume( const G4String& );
         G4double           GetWorldSizeX();
         G4double           GetWorldSizeYZ();
         G4double           GetCalorThickness();
         G4double           GetCalorSizeYZ();
         G4int              GetNbOfLayers();
   const G4VPhysicalVolume* GetphysiWorld();
   const G4VPhysicalVolume* GetAbsorber();
   const G4VPhysicalVolume* GetGap();
   

  private:
   SAXProcessor sxp;
   ProcessingConfigurator config;
   G4VPhysicalVolume* fWorld;
   ExN03CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
   ExN03DetectorMessenger* detectorMessenger;
};

#endif

