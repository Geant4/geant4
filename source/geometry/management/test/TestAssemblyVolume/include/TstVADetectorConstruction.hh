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
// $Id: TstVADetectorConstruction.hh,v 1.4 2001-07-11 09:59:24 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#ifndef TstVADetectorConstruction_h
#define TstVADetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class TstVADetectorMessenger;
class G4AssemblyVolume;

#include "g4std/vector"
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class TstVADetectorConstruction : public G4VUserDetectorConstruction
{
//  public:
    struct sClassic
    {
      G4LogicalVolume*                    caloLV;
      G4std::vector<G4VPhysicalVolume*>   PVs;
    };
    
  public:
    TstVADetectorConstruction();
    ~TstVADetectorConstruction();

  public:
    G4VPhysicalVolume*      Construct();
    void                    SwitchDetector();
    void                    SelectDetector(G4String val);
    void                    SelectMaterial(G4String val);

//  private:
    G4VPhysicalVolume*      SelectDetector();
    void                    SelectMaterialPointer();
    void                    ConstructClassic();
    void                    ConstructAssembly();
    void                    CleanClassic();
    void                    CleanAssembly();

  private:
    G4VPhysicalVolume*      worldVol;
    G4Material*             Air;
    G4Material*             Al;
    G4Material*             Pb;
    G4Material*             selectedMaterial;
    G4int                   detectorChoice;
    G4String                materialChoice;
    TstVADetectorMessenger* detectorMessenger;
     
  private:
    // Very private data
    G4LogicalVolume*        plateLV;
    sClassic                classicDetector;
    G4AssemblyVolume*       assemblyDetector;
};

#endif

