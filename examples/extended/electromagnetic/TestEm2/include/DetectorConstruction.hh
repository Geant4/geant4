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
// $Id: DetectorConstruction.hh,v 1.9 2006-06-29 16:49:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

class G4Tubs;
class G4LogicalVolume;
class G4UniformMagField;
class DetectorMessenger;

      const G4int MaxBin = 500;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:

     void SetMaterial(const G4String&);
     void SetLBining (G4ThreeVector);
     void SetRBining (G4ThreeVector);
     void SetMagField(G4double);

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();

     const
     G4VPhysicalVolume* GetEcal() {return physiEcal;};
     G4Material*    GetMaterial() {return myMaterial;};

     // Subdivision of absorber
     G4int    GetnLtot()          {return nLtot;};
     G4int    GetnRtot()          {return nRtot;};
     G4double GetdLradl()         {return dLradl;};
     G4double GetdRradl()         {return dRradl;};
     G4double GetdLlength()       {return dLlength;};
     G4double GetdRlength()       {return dRlength;};     
     G4double GetfullLength()     {return EcalLength;};
     G4double GetfullRadius()     {return EcalRadius;};

  private:

     void DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();

     G4int    nLtot,    nRtot;        // nb of bins: longitudinal and radial
     G4double dLradl,   dRradl;       // bin thickness (in radl unit)
     G4double dLlength, dRlength;     // bin thickness (in length unit)

     G4Material* myMaterial;          //pointer to the material
     G4UniformMagField* magField;     //pointer to the mag field

     G4double EcalLength;             //full length of the Calorimeter
     G4double EcalRadius;             //radius  of the Calorimeter

     G4Tubs*            solidEcal;    //pointer to the solid calorimeter
     G4LogicalVolume*   logicEcal;    //pointer to the logical calorimeter
     G4VPhysicalVolume* physiEcal;    //pointer to the physical calorimeter

     DetectorMessenger* detectorMessenger;  //pointer to the Messenger

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

