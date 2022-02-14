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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// and the DNA geometry given in the Geom_DNA example 
// shall cite the following Geant4-DNA collaboration publications:
// [1] NIM B 298 (2013) 47-54
// [2] Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

// Geant4
#include "globals.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "ChromosomeParameterisation.hh"

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::mm;
using CLHEP::degree;
using CLHEP::nanometer;
using CLHEP::micrometer;

static G4VisAttributes visInvBlue(false, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visInvWhite(false, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visInvPink(false, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visInvCyan(false, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visInvRed(false, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visInvGreen(false, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visWhite(true, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visPink(true, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visCyan(true, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visRed(true, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visGreen(true, G4Colour(0.0, 1.0, 0.0));


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() :
    G4VUserDetectorConstruction(),
    fBuildChromatineFiber(true),
    fBuildBases(false)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Water is defined from NIST material database
  G4NistManager *man = G4NistManager::Instance();
  man->FindOrBuildMaterial("G4_WATER");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::LoadChromosome(const char* filename,
                                          G4VPhysicalVolume* chromBox,
                                          G4LogicalVolume* logicBoxros)
{
  ChromosomeParameterisation* cp = new ChromosomeParameterisation(filename);
  new G4PVParameterised("box ros",
                        logicBoxros,
                        chromBox,
                        kUndefined,
                        cp->GetNumRosettes(),
                        cp);

  G4cout << filename << " done" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructDetector()
{
  if(fBuildBases == false && fBuildChromatineFiber == false)
  {
   G4cout <<"======================================================" << G4endl;
   G4cout <<"WARNING from DetectorConstruction::ConstructDetector:" << G4endl;
   G4cout << "As long as the flags fBuildBases and fBuildChromatineFiber are "
       "false, the output root file will be empty" << G4endl;
   G4cout << "This is intended for fast computation, display, testing ..."
          << G4endl;
   G4cout <<"======================================================" << G4endl;
  }

  G4String name;

  /***************************************************************************/
  //                               World
  /***************************************************************************/

  DefineMaterials();
  G4Material* waterMaterial = G4Material::GetMaterial("G4_WATER");

  G4Box* solidWorld = new G4Box("world", 10.0 * mm, 10.0 * mm, 10.0 * mm);
  G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,
                                                    waterMaterial,
                                                    "world");
  G4PVPlacement* physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);
  logicWorld->SetVisAttributes(&visInvWhite);

  /****************************************************************************/
  //                             Box nucleus
  /****************************************************************************/

  G4Box* solidTin = new G4Box("tin",
                              13 * micrometer,
                              10 * micrometer,
                              5 * micrometer);
  G4LogicalVolume* logicTin = new G4LogicalVolume(solidTin,
                                                  waterMaterial,
                                                  "tin");
  G4VPhysicalVolume* physiTin = new G4PVPlacement(0,
                                                  G4ThreeVector(),
                                                  "tin",
                                                  logicTin,
                                                  physiWorld,
                                                  false,
                                                  0);
  logicTin->SetVisAttributes(&visInvWhite);

  /****************************************************************************/
  //                            Cell nucleus
  /****************************************************************************/

  G4Ellipsoid* solidNucleus = new G4Ellipsoid("nucleus",
                                              11.82 * micrometer,
                                              8.52 * micrometer,
                                              3 * micrometer,
                                              0,
                                              0);
  G4LogicalVolume* logicNucleus = new G4LogicalVolume(solidNucleus,
                                                      waterMaterial,
                                                      "logic nucleus");
  G4VPhysicalVolume* physiNucleus = new G4PVPlacement(0,
                                   G4ThreeVector(),
                                   "physi nucleus",
                                   logicNucleus,
                                   physiTin,
                                   false,
                                   0);
  logicNucleus->SetVisAttributes(&visPink);

  /****************************************************************************/
  //                        Chromosomes territories
  /****************************************************************************/
  // NOTE: The only supported values for the rotation are
  // 0 and 90 degrees on the Y axis.
  G4double chromosomePositionSizeRotation[][7] = {
      {4.467, 2.835, 0, 1.557, 1.557, 1.557, 90},
      {-4.467, 2.835, 0, 1.557, 1.557, 1.557, 0},
      {4.423, -2.831, 0, 1.553, 1.553, 1.553, 90},
      {-4.423, -2.831, 0, 1.553, 1.553, 1.553, 0},
      {1.455, 5.63, 0, 1.455, 1.455, 1.455, 0},
      {-1.455, 5.63, 0, 1.455, 1.455, 1.455, 90},
      {1.435, 0, 1.392, 1.435, 1.435, 1.435, 0},
      {-1.435, 0, 1.392, 1.435, 1.435, 1.435, 90},
      {1.407, 0, -1.450, 1.407, 1.407, 1.407, 90}, // 5 right
      {-1.407, 0, -1.450, 1.407, 1.407, 1.407, 0}, // 5 left
      {1.380, -5.437, 0, 1.380, 1.380, 1.380, 0},
      {-1.380, -5.437, 0, 1.380, 1.380, 1.380, 90},
      {1.347, 2.782, -1.150, 1.347, 1.347, 1.347, 90},
      {-1.347, 2.782, -1.150, 1.347, 1.347, 1.347, 0},
      {1.311, -2.746, -1.220, 1.311, 1.311, 1.311, 90},
      {-1.311, -2.746, -1.220, 1.311, 1.311, 1.311, 0},
      {7.251, -2.541, 0, 1.275, 1.275, 1.275, 0},
      {-6.701, 0, -0.85, 1.275, 1.275, 1.275, 90},
      {4.148, 0, 1.278, 1.278, 1.278, 1.278, 90}, // 10 right
      {-4.148, 0, 1.278, 1.278, 1.278, 1.278, 0}, // 10 left
      {4.147, 0, -1.277, 1.277, 1.277, 1.277, 0},
      {-4.147, 0, -1.277, 1.277, 1.277, 1.277, 90},
      {8.930, 0.006, 0, 1.272, 1.272, 1.272, 90},
      {-7.296, 2.547, 0, 1.272, 1.272, 1.272, 90},
      {1.207, -2.642, 1.298, 1.207, 1.207, 1.207, 0},
      {-1.207, -2.642, 1.298, 1.207, 1.207, 1.207, 90},
      {1.176, 2.611, 1.368, 1.176, 1.176, 1.176, 0},
      {-1.176, 2.611, 1.368, 1.176, 1.176, 1.176, 90},
      {4.065, 5.547, 0, 1.155, 1.155, 1.155, 90}, // 15 right
      {-4.065, 5.547, 0, 1.155, 1.155, 1.155, 0}, // 15 left
      {6.542, 0.159, 1.116, 1.116, 1.116, 1.116, 0},
      {-9.092, 0, 0, 1.116, 1.116, 1.116, 0},
      {6.507, 0.159, -1.081, 1.081, 1.081, 1.081, 90},
      {-7.057, -2.356, 0, 1.081, 1.081, 1.081, 90},
      {3.824, -5.448, 0, 1.064, 1.064, 1.064, 90},
      {-3.824, -5.448, 0, 1.064, 1.064, 1.064, 0},
      {5.883, -5.379, 0, 0.995, 0.995, 0.995, 0},
      {-9.133, -2.111, 0, 0.995, 0.995, 0.995, 0},
      {6.215, 5.387, 0, 0.995, 0.995, 0.995, 0}, // 20 right
      {-6.971, -4.432, 0, 0.995, 0.995, 0.995, 90}, // 20 left
      {9.583, 2.177, 0, 0.899, 0.899, 0.899, 90},
      {-9.467, 2.03, 0, 0.899, 0.899, 0.899, 0},
      {9.440, -2.180, 0, 0.914, 0.914, 0.914, 90},
      {-6.34, 0, 1.339, 0.914, 0.914, 0.914, 0},
      {-6.947, 4.742, 0, 0.923, 0.923, 0.923, 90}, // Y
      {7.354, 2.605, 0, 1.330, 1.330, 1.330, 0} // X
  };

  G4RotationMatrix* rotch = new G4RotationMatrix;
  rotch->rotateY(90 * degree);

  vector<G4VPhysicalVolume*> physiBox(48);

  for (unsigned int i = 0; i < countof(chromosomePositionSizeRotation); i++)
  {
    G4double* p = &chromosomePositionSizeRotation[i][0];
    G4double* size = &chromosomePositionSizeRotation[i][3];
    G4double rotation = chromosomePositionSizeRotation[i][6];
    G4ThreeVector pos(p[0] * micrometer, p[1] * micrometer, p[2] * micrometer);
    G4RotationMatrix* rot = rotation == 0 ? 0 : rotch;

    ostringstream ss;
    ss << "box" << (i / 2) + 1 << (i % 2 ? 'l' : 'r');
    name = ss.str();
    ss.str("");
    ss.clear();

    /*
     snprintf(name, countof(name), "box%d%c",
     (i / 2) + 1, i % 2 ? 'l' : 'r');
     */
    G4Box* solidBox = new G4Box(name,
                                size[0] * micrometer,
                                size[1] * micrometer,
                                size[2] * micrometer);
    G4LogicalVolume* logicBox = new G4LogicalVolume(solidBox,
                                                    waterMaterial,
                                                    name);
    physiBox[i] = new G4PVPlacement(rot,
                                    pos,
                                    "chromo",
                                    logicBox,
                                    physiNucleus,
                                    false,
                                    0);
    logicBox->SetVisAttributes(&visBlue);
  }

  /**************************************************************************/
  //                 Box containing the chromatin flowers
  /**************************************************************************/

  G4Tubs* solidBoxros = new G4Tubs("solid box ros",
                                   0 * nanometer,
                                   399 * nanometer,
                                   20 * nanometer,
                                   0 * degree,
                                   360 * degree);
  G4LogicalVolume* logicBoxros = new G4LogicalVolume(solidBoxros,
                                                     waterMaterial,
                                                     "box ros");
  logicBoxros->SetVisAttributes(&visInvBlue);

  //Loading flower box position for each chromosome territory

  for (int k = 0; k < 22; k++)
  {
    ostringstream oss;
    oss << "chromo" << k + 1 << ".dat";
    name = oss.str();
    oss.str("");
    oss.clear();
    //snprintf(name, countof(name), "chromo%d.dat", k + 1);
    LoadChromosome(name.c_str(), physiBox[k * 2], logicBoxros);
    LoadChromosome(name.c_str(), physiBox[k * 2 + 1], logicBoxros);
  }

  LoadChromosome("chromoY.dat", physiBox[44], logicBoxros);
  LoadChromosome("chromoX.dat", physiBox[45], logicBoxros);

  /****************************************************************************/
  if (fBuildChromatineFiber)
  {
    // chromatin fiber envelope
    G4Tubs* solidEnv = new G4Tubs("chromatin fiber",
                                  0,
                                  15.4 * nanometer,
                                  80.5 * nanometer,
                                  0 * degree,
                                  360 * degree);
    G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv,
                                                    waterMaterial,
                                                    "LV chromatin fiber");
    logicEnv->SetVisAttributes(&visInvPink);

    // Chromatin fiber position
    for (G4int i = 0; i < 7; i++)
    {
      G4RotationMatrix* rotFiber = new G4RotationMatrix;
      rotFiber->rotateX(90 * degree);
      rotFiber->rotateY(i * 25.72 * degree);
      G4ThreeVector posFiber = G4ThreeVector(0, 152 * nanometer, 0);
      posFiber.rotateZ(i * 25.72 * degree);
      new G4PVPlacement(rotFiber,
                        posFiber,
                        logicEnv,
                        "physi env",
                        logicBoxros,
                        false,
                        0);

      rotFiber = new G4RotationMatrix;
      rotFiber->rotateX(90 * degree);
      rotFiber->rotateY((7 + i) * 25.72 * degree);
      posFiber = G4ThreeVector(0, 152 * nanometer, 0);
      posFiber.rotateZ((7 + i) * 25.72 * degree);
      new G4PVPlacement(rotFiber,
                        posFiber,
                        logicEnv,
                        "physi env",
                        logicBoxros,
                        false,
                        0);

      rotFiber = new G4RotationMatrix;
      rotFiber->rotateX(90 * degree);
      rotFiber->rotateY((25.72 + (i - 14) * 51.43) * degree);
      posFiber = G4ThreeVector(-36.5 * nanometer, 312 * nanometer, 0);
      posFiber.rotateZ((i - 14) * 51.43 * degree);
      new G4PVPlacement(rotFiber,
                        posFiber,
                        logicEnv,
                        "physi env",
                        logicBoxros,
                        false,
                        0);

      rotFiber = new G4RotationMatrix;
      rotFiber->rotateX(90 * degree);
      rotFiber->rotateY(180 * degree);
      rotFiber->rotateY((i - 21) * 51.43 * degree);
      posFiber = G4ThreeVector(-103 * nanometer, 297 * nanometer, 0);
      posFiber.rotateZ((i - 21) * 51.43 * degree);
      new G4PVPlacement(rotFiber,
                        posFiber,
                        logicEnv,
                        "physi env",
                        logicBoxros,
                        false,
                        0);

    }

    if (fBuildBases)
    {
      // Histones
      G4Tubs* solidHistone = new G4Tubs("solid histone",
                                        0,
                                        3.25 * nanometer,
                                        2.85 * nanometer,
                                        0 * degree,
                                        360 * degree);
      G4LogicalVolume* logicHistone = new G4LogicalVolume(solidHistone,
                                                          waterMaterial,
                                                          "logic histone");

      //Base pair
      G4Orb* solidBp1 = new G4Orb("blue sphere", 0.17 * nanometer);
      G4LogicalVolume* logicBp1 = new G4LogicalVolume(solidBp1,
                                                      waterMaterial,
                                                      "logic blue sphere");
      G4Orb* solidBp2 = new G4Orb("pink sphere", 0.17 * nanometer);
      G4LogicalVolume* logicBp2 = new G4LogicalVolume(solidBp2,
                                                      waterMaterial,
                                                      "logic pink sphere");

      //Phosphodiester group

      G4Orb* solidSugar_48em1_nm = new G4Orb("sugar", 0.48 * nanometer);

      G4ThreeVector posi(0.180248 * nanometer,
                         0.32422 * nanometer,
                         0.00784 * nanometer);
      G4UnionSolid* uniDNA = new G4UnionSolid("move",
                                              solidSugar_48em1_nm,
                                              solidSugar_48em1_nm,
                                              0,
                                              posi);

      G4ThreeVector posi2(-0.128248 * nanometer,
                          0.41227 * nanometer,
                          0.03584 * nanometer);
      G4UnionSolid* uniDNA2 = new G4UnionSolid("move2",
                                               solidSugar_48em1_nm,
                                               solidSugar_48em1_nm,
                                               0,
                                               posi2);

      /************************************************************************
       Phosphodiester group Position
       ************************************************************************/

      for (G4int n = 2; n < 200; n++)
      {
      G4double SP1[2][3] = {
          { (-0.6 * nanometer) * cos(n * 0.26),
            0, (0.6* nanometer) * sin(n * 0.26) },
          { (0.6 * nanometer) * cos(n * 0.26),
            0, (-0.6 * nanometer) * sin(0.26 * n) }
      };
      G4double matriceSP1[3][3] = {
          { cos(n * 0.076), -sin(n * 0.076), 0 },
          { sin(n * 0.076), cos(n * 0.076), 0 },
          { 0, 0, 1 }
      };
      G4double matriceSP2[2][3];

      for (G4int i = 0; i < 3; i++)
      {
        G4double sumSP1 = 0;
        G4double sumSP2 = 0;
        for (G4int j = 0; j < 3; j++)
        {
          sumSP1 += matriceSP1[i][j] * SP1[0][j];
          sumSP2 += matriceSP1[i][j] * SP1[1][j];
        }
        matriceSP2[0][i] = sumSP1;
        matriceSP2[1][i] = sumSP2;
      }

      G4double heliceSP[3] = {
          (4.85 * nanometer) * cos(n * 0.076),
          (4.85 * nanometer) * sin(n * 0.076),
          (n * 0.026 * nanometer)
      };

      for (G4int i = 0; i < 3; i++)
      {
        matriceSP2[0][i] += heliceSP[i];
        matriceSP2[1][i] += heliceSP[i];
      }
      G4ThreeVector posSugar1(matriceSP2[0][2],
                              matriceSP2[0][1],
                              (matriceSP2[0][0]) - (4.25 * nanometer));
      G4ThreeVector posSugar2(matriceSP2[1][2],
                              matriceSP2[1][1],
                              (matriceSP2[1][0]) - (5.45 * nanometer));

      ostringstream ss;
      ss << "sugar_" << n;
      name = ss.str().c_str();
      ss.str("");
      ss.clear();

      //  snprintf(name, countof(name), "sugar %d", n);
      uniDNA = new G4UnionSolid(name,
                                uniDNA,
                                solidSugar_48em1_nm,
                                0,
                                posSugar1);

      ss << "sugar_" << n;
      name = ss.str().c_str();
      ss.str("");
      ss.clear();

      //  snprintf(name, countof(name), "sugar %d", n);
      uniDNA2 = new G4UnionSolid(name,
                                 uniDNA2,
                                 solidSugar_48em1_nm,
                                 0,
                                 posSugar2);
    }
    G4LogicalVolume* logicSphere3 = new G4LogicalVolume(uniDNA,
                                                        waterMaterial,
                                                        "logic sugar 2");
    G4LogicalVolume* logicSphere4 = new G4LogicalVolume(uniDNA2,
                                                        waterMaterial,
                                                        "logic sugar 4");

    /**************************************************************************
     Base pair Position
     **************************************************************************/
      for (G4int n = 0; n < 200; n++)
      {
        G4double bp1[2][3] = {
            { (-0.34 * nanometer) * cos(n * 0.26),
            0, (0.34* nanometer) * sin(n * 0.26) },
            { (0.34 * nanometer) * cos(n * 0.26),
            0, (-0.34 * nanometer) * sin(0.26 * n) }
        };
        G4double matriceBP1[3][3] = {
            { cos(n * 0.076), -sin(n * 0.076), 0 },
            {sin(n * 0.076), cos(n * 0.076), 0 },
            { 0, 0, 1 }
        };
        G4double matriceBP2[2][3];

        for (G4int i = 0; i < 3; i++)
        {
          G4double sumBP1 = 0;
          G4double sumBP2 = 0;
          for (G4int j = 0; j < 3; j++)
          {
            sumBP1 += matriceBP1[i][j] * bp1[0][j];
            sumBP2 += matriceBP1[i][j] * bp1[1][j];
          }
          matriceBP2[0][i] = sumBP1;
          matriceBP2[1][i] = sumBP2;
        }
        G4double heliceBP[3] = {
            (4.8 * nanometer) * cos(n * 0.076),
            (4.8 * nanometer) * sin(n * 0.076),
            n * 0.026 * nanometer
        };

        for (G4int i = 0; i < 3; i++)
        {
          matriceBP2[0][i] += heliceBP[i];
          matriceBP2[1][i] += heliceBP[i];
        }
        G4ThreeVector position1(matriceBP2[0][2],
                                matriceBP2[0][1],
                                matriceBP2[0][0] - (4.25 * nanometer));
        G4ThreeVector position2(matriceBP2[1][2],
                                matriceBP2[1][1],
                                matriceBP2[1][0] - (5.45 * nanometer));

        new G4PVPlacement(0,
                          position1,
                          logicBp1,
                          "physi blue sphere",
                          logicSphere3,
                          false,
                          0);
        new G4PVPlacement(0,
                          position2,
                          logicBp2,
                          "physi pink sphere",
                          logicSphere4,
                          false,
                          0);
        }

  /****************************************************************************/
  //                 Initial position of different elements
  /****************************************************************************/
      // DNA and histone positions
      for (int j = 0; j < 90; j++)
      {
        // DNA (bp-SP)
        G4RotationMatrix* rotStrand1 = new G4RotationMatrix;
        rotStrand1->rotateZ(j * -51.43 * degree);
        G4ThreeVector posStrand1(-2.7 * nanometer,
                                 9.35 * nanometer,
                                 (-69.9 * nanometer) + (j * 1.67 * nanometer));
        posStrand1.rotateZ(j * 51.43 * degree);
        new G4PVPlacement(rotStrand1,
                          posStrand1,
                          logicSphere3,
                          "physi sugar 2",
                          logicEnv,
                          false,
                          0);

        G4RotationMatrix* rotStrand2 = new G4RotationMatrix;
        rotStrand2->rotateZ(j * -51.43 * degree);
        G4ThreeVector posStrand2(-2.7 * nanometer,
                                 9.35 * nanometer,
                                 (-68.7 * nanometer) + (j * 1.67 * nanometer));
        posStrand2.rotateZ(j * 51.43 * degree);
        new G4PVPlacement(rotStrand2,
                          posStrand2,
                          logicSphere4,
                          "physi sugar 4",
                          logicEnv,
                          false,
                          0);

        // histones
        G4RotationMatrix* rotHistone = new G4RotationMatrix;
        rotHistone->rotateY(90 * degree);
        rotHistone->rotateX(j * (-51.43 * degree));
        G4ThreeVector posHistone(0.0,
                                 9.35 * nanometer,
                                 (-74.15 + j * 1.67) * nanometer);
        posHistone.rotateZ(j * 51.43 * degree);
        new G4PVPlacement(rotHistone,
                          posHistone,
                          logicHistone,
                          "PV histone",
                          logicEnv,
                          false,
                          0);
      }
      /************************************************************************/
      //                        Visualisation colors
      /************************************************************************/

      logicBp1->SetVisAttributes(&visInvCyan);
      logicBp2->SetVisAttributes(&visInvPink);

      logicSphere3->SetVisAttributes(&visInvWhite);
      logicSphere4->SetVisAttributes(&visInvRed);

      logicHistone->SetVisAttributes(&visInvBlue);
    }
  }

  G4cout << "Geometry has been loaded" << G4endl;
  return physiWorld;
}
