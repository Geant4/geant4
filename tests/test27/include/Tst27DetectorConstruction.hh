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
// $Id: Tst27DetectorConstruction.hh,v 1.3 2004-04-07 15:54:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst27DetectorConstruction_h
#define Tst27DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst27DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst27DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst27DetectorConstruction();
    ~Tst27DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SelectMaterial(G4String val);

  private:
     void SelectMaterialPointer();

     G4LogicalVolume*   simpleBoxLog;
     G4Material* theH;
     G4Material* theSi;
     G4Material* theCu;
     G4Material* thePb;
     G4Material* theU;
     G4Material* theH2O;
     
     G4Material* selectedMaterial;
     G4String materialChoice;
     Tst27DetectorMessenger * detectorMessenger;
};

#endif

