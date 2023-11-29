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
// Authors: J. Naoki D. Kondo (UCSF, US) : 10/10/2021 
//          J. Ramos-Mendez and B. Faddegon (UCSF, US)
//
/// \file PhysGeoImport.cc
/// \brief Implementation of the plasmid load methods for the geometry

#include "PhysGeoImport.hh"
#include "G4DNAChemistryManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Ellipsoid.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

PhysGeoImport::PhysGeoImport() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LogicalVolume* PhysGeoImport::CreateLogicVolumeXYZ(G4String fileName) {
  G4NistManager * man = G4NistManager::Instance();
  fEnvelopeWater = man->FindOrBuildMaterial("G4_WATER");
  fpWater = man->BuildMaterialWithNewDensity("G4_WATER_MODIFIED",
                                       "G4_WATER",1.0*g/cm/cm/cm);

  ReadFile(fileName);

  G4double des1 = 1.1344640137963142;
  G4double des2 = des1 + (CLHEP::pi*.5);
  G4double ang = 0.6283185307179586;
  G4double bet1 = 0.6283185307179586 * 2;
  G4double posi = 1.0471975511965976;
  G4double sep = .1*angstrom;
  //Geometries Sizes
  G4double PxRs = 2.9389169420478556 * angstrom; 
  G4double PyRs = 2.9389169420478556 * angstrom;
  G4double PzRs = 2.9389169420478556 * angstrom;
  G4double PxYs = 2.7 * angstrom;
  G4double PyYs = 2.7 * angstrom;
  G4double PzYs = 2.7 * angstrom;
  G4double PxBp = 2.45 * angstrom;
  G4double PyBp = 2.45 * angstrom;
  G4double PzBp = 2.45 * angstrom;
 

  G4double xin = -170 * angstrom;
  G4double yin = -170 * angstrom;
  G4double zin = -170 * angstrom;
  G4double xfn =  170 * angstrom;
  G4double yfn =  170 * angstrom;
  G4double zfn =  170 * angstrom;  

  G4int nVertex = fVertexes.size();

  // Envelope

  std::string boxNameSolid = fGeoName + "_solid";
  G4Box* box_solid = new G4Box(boxNameSolid, 0.5*(fXMax-fXMin)+0.5*3.4*nm,
                               0.5*(fYMax-fYMin)+0.5*3.4*nm,
                               0.5*(fZMax-fZMin)+0.5*3.4*nm);

  G4String boxNameLogic       = fGeoName + "_logic";
  G4LogicalVolume* box_logic  = new G4LogicalVolume(box_solid, 
                                                    fEnvelopeWater, 
                                                    boxNameLogic, 0,0,0);

  //Desoxyribose

  G4Ellipsoid *RSolidSugar = new G4Ellipsoid("sdeoxyribose", 
                                             PxRs, 
                                             PyRs, 
                                             PzRs, -PzRs, .445*PzRs);
  G4LogicalVolume* RSugar  = new G4LogicalVolume(RSolidSugar, 
                                                 fpWater,
                                                 "ldeoxyribose",
                                                 0,0,0);
  G4VisAttributes* MyVisAtt_Rs = new G4VisAttributes(G4Colour(G4Colour::Red()));
  MyVisAtt_Rs->SetForceSolid(true);
  RSugar->SetVisAttributes(MyVisAtt_Rs);

  //Phosphoric Acid

  G4Ellipsoid *YSolidSugar = new G4Ellipsoid("sphosphate", 
                                             PxYs, 
                                             PyYs, 
                                             PzYs, 
                                             -PzYs, 
                                             .9 * angstrom);

  G4LogicalVolume* YSugar = new G4LogicalVolume(YSolidSugar, 
                                                fpWater, 
                                                "lphosphate",
                                                0,0,0);
  G4VisAttributes* MyVisAtt_Ys = new G4VisAttributes(G4Colour(G4Colour::Yellow()));
  MyVisAtt_Ys->SetForceSolid(true);
  YSugar->SetVisAttributes(MyVisAtt_Ys);

  //Base Pairs

  G4Ellipsoid* Base1a = new G4Ellipsoid("BasePair1a", 
                                        PxBp, 
                                        PyBp, 
                                        PzBp, 
                                        -PzBp, 
                                        1.15 * angstrom);
  G4LogicalVolume* BaseP1a = new G4LogicalVolume(Base1a, fpWater, "BasePair1a", 0,0,0);
  G4VisAttributes* MyVisAtt_Bp1a = new G4VisAttributes(G4Colour(G4Colour::Green()));
  MyVisAtt_Bp1a->SetForceSolid(true);
  BaseP1a->SetVisAttributes(MyVisAtt_Bp1a);

  G4Ellipsoid* Base1b= new G4Ellipsoid("BasePair1b", 
                                       PxBp, 
                                       PyBp, 
                                       PzBp, 
                                       1.15 * angstrom, 
                                       PzBp);
  G4LogicalVolume* BaseP1b = new G4LogicalVolume(Base1b, fpWater, "BasePair1b", 0,0,0);

  G4VisAttributes* MyVisAtt_Bp1b = new G4VisAttributes(G4Colour(G4Colour::Green()));
  MyVisAtt_Bp1b->SetForceSolid(true);
  BaseP1b->SetVisAttributes(MyVisAtt_Bp1b);

  G4int index = 0; 
  G4double cAngle = 0;
  G4double pi =  CLHEP::pi;
  for (int vertex = 0; vertex < nVertex - 1; vertex++ ) {
    xin = fVertexes[vertex][0]-fOffsetX;
    yin = fVertexes[vertex][1]-fOffsetY;
    zin = fVertexes[vertex][2]-fOffsetZ;
    xfn = fVertexes[vertex+1][0]-fOffsetX;
    yfn = fVertexes[vertex+1][1]-fOffsetY;
    zfn = fVertexes[vertex+1][2]-fOffsetZ;

    G4double phi0   = std::atan2(zfn-zin,std::sqrt(((xfn-xin)*(xfn-xin))+
                                                       ((yfn-yin)*(yfn-yin))));
    G4double theta0 = std::atan2(xfn-xin,yfn-yin);
    G4double lenght = std::sqrt(((xfn-xin)*(xfn-xin))+
                                    ((yfn-yin)*(yfn-yin))+
                                    ((zfn-zin)*(zfn-zin)));
    G4double dl = 1.0 / (lenght / (3.4*angstrom));
    G4int nChain = (fVertexes[vertex] - fVertexes[vertex+1]).mag() / (0.34 * nm);

    for (G4int nseg = 0; nseg < nChain ; nseg++) {
      cAngle += ang;

      G4double theta = cAngle;
      G4double x1 = 0;
      G4double y1 = 0;
      G4double z1 = ((2*PzBp) + PzRs + sep);

      G4double x2 = 0; 
      G4double y2 = 0; 
      G4double z2 = ((2*PzBp) + PzRs + sep);

      G4ThreeVector plus2 = G4ThreeVector(0,0,(.5 * PzRs) + PzYs);
      plus2.rotateX(-des1);
      plus2.rotateZ(-posi);

      G4ThreeVector plus2alt = G4ThreeVector(0,0,(.5 * PzRs) + PzYs);
      plus2alt.rotateX(-des1);
      plus2alt.rotateZ(posi);

      G4double x3 = 0;
      G4double y3 = 0;
      G4double z3 = PzBp + sep;

      G4ThreeVector position1i = G4ThreeVector(x1,y1,z1);
      G4ThreeVector position2i = G4ThreeVector(x2,y2,z2)+plus2;
      G4ThreeVector position2ialt = G4ThreeVector(x2,y2,-z2)-plus2alt;
      G4ThreeVector position3i = G4ThreeVector(x3,y3,z3);

      position1i.rotateY(theta);
      position2i.rotateY(theta);
      position2ialt.rotateY(theta);
      position3i.rotateY(theta);

      G4double x = dl*nseg*(xfn-xin)+xin;
      G4double y = dl*nseg*(yfn-yin)+yin;
      G4double z = dl*nseg*(zfn-zin)+zin;

      position1i.rotateX(phi0);
      position2i.rotateX(phi0);
      position2ialt.rotateX(phi0);
      position3i.rotateX(phi0);

      position1i.rotateZ(-theta0);
      position2i.rotateZ(-theta0);
      position2ialt.rotateZ(-theta0);
      position3i.rotateZ(-theta0);

      G4double yrot1 = theta;
      G4double xrot1 = -des1;
      G4RotationMatrix rotm1 = G4RotationMatrix();
      rotm1.rotateX(xrot1);
      rotm1.rotateZ(-posi);
      rotm1.rotateY(yrot1);
      rotm1.rotateX(phi0);
      rotm1.rotateZ(-theta0);
      G4ThreeVector position1 = position1i + G4ThreeVector(x,y,z);
      G4Transform3D transform1(rotm1, position1);
  
      G4double yrot1alt = theta + pi;
      G4double xrot1alt = des1;
      G4RotationMatrix rotm1alt = G4RotationMatrix();
      rotm1alt.rotateX(xrot1alt);
      rotm1alt.rotateZ(-posi);
      rotm1alt.rotateY(yrot1alt);
      rotm1alt.rotateX(phi0);
      rotm1alt.rotateZ(-theta0);
      G4ThreeVector position1alt = -position1i + G4ThreeVector(x,y,z);
      G4Transform3D transform1alt(rotm1alt,position1alt);

      G4double yrot2 = theta;
      G4double xrot2 = -des2;
      G4RotationMatrix rotm2 = G4RotationMatrix();
      rotm2.rotateX(xrot2);
      rotm2.rotateY(yrot2 - bet1 + 0.8726646259971648);
      rotm2.rotateX(phi0);
      rotm2.rotateZ(-theta0);
      G4ThreeVector position2 = position2i + G4ThreeVector(x,y,z);
      G4Transform3D transform2(rotm2,position2);

      G4double yrot2alt = theta + pi;
      G4double xrot2alt = des2;
      G4RotationMatrix rotm2alt = G4RotationMatrix();
      rotm2alt.rotateX(xrot2alt);
      rotm2alt.rotateY(yrot2alt + bet1 - 0.8726646259971648);
      rotm2alt.rotateX(phi0);
      rotm2alt.rotateZ(-theta0);
      G4ThreeVector position2alt = position2ialt + G4ThreeVector(x,y,z);
      G4Transform3D transform2alt(rotm2alt,position2alt);

      G4double yrot3 = theta;
      G4RotationMatrix rotm3 = G4RotationMatrix();
      rotm3.rotateX(-pi/2);
      rotm3.rotateZ(-ang);
      rotm3.rotateY(yrot3);
      rotm3.rotateX(phi0);
      rotm3.rotateZ(-theta0);
      G4ThreeVector position3 = position3i + G4ThreeVector(x,y,z);
      G4Transform3D transform3(rotm3,position3);

      G4double yrot3alt = theta + pi;
      G4RotationMatrix rotm3alt = G4RotationMatrix();
      rotm3alt.rotateX(pi/2);
      rotm3alt.rotateZ(-ang);
      rotm3alt.rotateY(yrot3alt);
      rotm3alt.rotateX(phi0);
      rotm3alt.rotateZ(-theta0);
      G4ThreeVector position3alt = -position3i + G4ThreeVector(x,y,z);
      G4Transform3D transform3alt(rotm3alt,position3alt);

     
      new G4PVPlacement(transform1,
                        RSugar,
                        "deoxyribose1",
                        box_logic,
                        false,
                        index,
                        false);

      new G4PVPlacement(transform1alt,
                        RSugar,
                        "deoxyribose2",
                        box_logic,
                        false,
                        index,
                        false);

      new G4PVPlacement(transform3,
                        BaseP1a,
                        "BasePair1",
                        box_logic,
                        false,
                        index,
                        false);

      new G4PVPlacement(transform3alt,
                        BaseP1a,
                        "BasePair2",
                        box_logic,
                        false,
                        index,
                        false);

      new G4PVPlacement(transform2,
                        YSugar,
                        "phosphate1",
                        box_logic,
                        false,
                        index,
                        false);

      new G4PVPlacement(transform2alt,
                        YSugar,
                        "phosphate2",
                        box_logic,
                        false,
                        index,
                        false);

      G4ThreeVector Deoxy1 = position1;
      G4ThreeVector Deoxy2 = position1alt;

      fSampleDNAPositions.push_back(Deoxy1);
      fSampleDNAPositions.push_back(Deoxy2);
      fSampleDNANames.push_back("Deoxyribose");
      fSampleDNANames.push_back("Deoxyribose");
      fSampleDNADetails.push_back({-1,index,1});
      fSampleDNADetails.push_back({-1,index,2});

      index++;
    
    }
  }

  return box_logic;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void PhysGeoImport::ReadFile(G4String fileName) {
  G4double x, y, z;
  fXMin = 1*mm, fYMin = 1*mm, fZMin = 1*mm;
  fXMax =-1*mm, fYMax =-1*mm, fZMax =-1*mm;
  std::ifstream plasmidFile(fileName);
  while(true) {
    plasmidFile >> x >> y >> z;
    if ( !plasmidFile.good() ) break;
    x *= nm;
    y *= nm;
    z *= nm;
    fVertexes.push_back(G4ThreeVector(x, y, z));
    if ( fXMin > x ) fXMin = x;
    if ( fXMax < x ) fXMax = x;
    if ( fYMin > y ) fYMin = y;
    if ( fYMax < y ) fYMax = y;
    if ( fZMin > z ) fZMin = z;
    if ( fZMax < z ) fZMax = z;
  }
  plasmidFile.close();
  fOffsetX = (fXMin + fXMax)*0.5;
  fOffsetY = (fYMin + fYMax)*0.5;
  fOffsetZ = (fZMin + fZMax)*0.5;

  std::vector<G4ThreeVector> VertRed;
  for(size_t i = 0; i < fVertexes.size(); i++) {
    if (i % 15 == 0)
      VertRed.push_back(fVertexes[i]);
  }

  VertRed.push_back(fVertexes[0]);

  fVertexes = VertRed;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
