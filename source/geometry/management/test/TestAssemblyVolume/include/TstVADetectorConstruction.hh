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
// $Id: TstVADetectorConstruction.hh,v 1.6 2006-06-29 18:34:15 gunter Exp $
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

#include <vector>
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class TstVADetectorConstruction : public G4VUserDetectorConstruction
{
//  public:
    struct sClassic
    {
      G4LogicalVolume*                    caloLV;
      std::vector<G4VPhysicalVolume*>   PVs;
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

