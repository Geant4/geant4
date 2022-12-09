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
//

#include "G4ScoringBox.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4VScoreColorMap.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4Polyhedron.hh"

#include "G4ScoringManager.hh"
#include "G4StatDouble.hh"

#include "G4SystemOfUnits.hh"

#include <map>
#include <fstream>

G4ScoringBox::G4ScoringBox(G4String wName)
  : G4VScoringMesh(wName)
  , fSegmentDirection(-1)
{
  fShape                = MeshShape::box;
  fDivisionAxisNames[0] = "X";
  fDivisionAxisNames[1] = "Y";
  fDivisionAxisNames[2] = "Z";
}

void G4ScoringBox::SetupGeometry(G4VPhysicalVolume* fWorldPhys)
{
  if(verboseLevel > 9)
    G4cout << "G4ScoringBox::SetupGeometry() ..." << G4endl;

  // World
  G4VPhysicalVolume* scoringWorld = fWorldPhys;
  G4LogicalVolume* worldLogical   = scoringWorld->GetLogicalVolume();

  // Scoring Mesh
  if(verboseLevel > 9)
    G4cout << fWorldName << G4endl;
  G4String boxName = fWorldName;

  if(verboseLevel > 9)
    G4cout << fSize[0] << ", " << fSize[1] << ", " << fSize[2] << G4endl;
  G4VSolid* boxSolid = new G4Box(boxName + "0", fSize[0], fSize[1], fSize[2]);
  auto  boxLogical =
    new G4LogicalVolume(boxSolid, nullptr, boxName + "_0");
  new G4PVPlacement(fRotationMatrix, fCenterPosition, boxLogical, boxName + "0",
                    worldLogical, false, 0);

  G4String layerName[2] = { boxName + "_1", boxName + "_2" };
  G4VSolid* layerSolid[2];
  G4LogicalVolume* layerLogical[2];

  //-- fisrt nested layer (replicated to x direction)
  if(verboseLevel > 9)
    G4cout << "layer 1 :" << G4endl;
  layerSolid[0] =
    new G4Box(layerName[0], fSize[0] / fNSegment[0], fSize[1], fSize[2]);
  layerLogical[0] = new G4LogicalVolume(layerSolid[0], nullptr, layerName[0]);
  if(fNSegment[0] > 1)
  {
    if(verboseLevel > 9)
      G4cout << "G4ScoringBox::Construct() : Replicate to x direction"
             << G4endl;
    if(G4ScoringManager::GetReplicaLevel() > 0)
    {
      new G4PVReplica(layerName[0], layerLogical[0], boxLogical, kXAxis,
                      fNSegment[0], fSize[0] / fNSegment[0] * 2.);
    }
    else
    {
      new G4PVDivision(layerName[0], layerLogical[0], boxLogical, kXAxis,
                       fNSegment[0], 0.);
    }
  }
  else if(fNSegment[0] == 1)
  {
    if(verboseLevel > 9)
      G4cout << "G4ScoringBox::Construct() : Placement" << G4endl;
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), layerLogical[0],
                      layerName[0], boxLogical, false, 0);
  }
  else
    G4cerr << "ERROR : G4ScoringBox::SetupGeometry() : invalid parameter ("
           << fNSegment[0] << ") "
           << "in placement of the first nested layer." << G4endl;

  if(verboseLevel > 9)
  {
    G4cout << fSize[0] / fNSegment[0] << ", " << fSize[1] << ", " << fSize[2]
           << G4endl;
    G4cout << layerName[0] << ": kXAxis, " << fNSegment[0] << ", "
           << 2. * fSize[0] / fNSegment[0] << G4endl;
  }

  // second nested layer (replicated to y direction)
  if(verboseLevel > 9)
    G4cout << "layer 2 :" << G4endl;
  layerSolid[1]   = new G4Box(layerName[1], fSize[0] / fNSegment[0],
                            fSize[1] / fNSegment[1], fSize[2]);
  layerLogical[1] = new G4LogicalVolume(layerSolid[1], nullptr, layerName[1]);
  if(fNSegment[1] > 1)
  {
    if(verboseLevel > 9)
      G4cout << "G4ScoringBox::Construct() : Replicate to y direction"
             << G4endl;
    if(G4ScoringManager::GetReplicaLevel() > 1)
    {
      new G4PVReplica(layerName[1], layerLogical[1], layerLogical[0], kYAxis,
                      fNSegment[1], fSize[1] / fNSegment[1] * 2.);
    }
    else
    {
      new G4PVDivision(layerName[1], layerLogical[1], layerLogical[0], kYAxis,
                       fNSegment[1], 0.);
    }
  }
  else if(fNSegment[1] == 1)
  {
    if(verboseLevel > 9)
      G4cout << "G4ScoringBox::Construct() : Placement" << G4endl;
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), layerLogical[1],
                      layerName[1], layerLogical[0], false, 0);
  }
  else
    G4cerr << "ERROR : G4ScoringBox::SetupGeometry() : invalid parameter ("
           << fNSegment[1] << ") "
           << "in placement of the second nested layer." << G4endl;

  if(verboseLevel > 9)
  {
    G4cout << fSize[0] / fNSegment[0] << ", " << fSize[1] / fNSegment[1] << ", "
           << fSize[2] << G4endl;
    G4cout << layerName[1] << ": kYAxis, " << fNSegment[1] << ", "
           << 2. * fSize[1] / fNSegment[1] << G4endl;
  }

  // mesh elements (replicated to z direction)
  if(verboseLevel > 9)
    G4cout << "mesh elements :" << G4endl;
  G4String elementName = boxName + "_3";
  G4VSolid* elementSolid =
    new G4Box(elementName, fSize[0] / fNSegment[0], fSize[1] / fNSegment[1],
              fSize[2] / fNSegment[2]);
  fMeshElementLogical = new G4LogicalVolume(elementSolid, nullptr, elementName);
  if(fNSegment[2] > 1)
  {
    if(verboseLevel > 9)
      G4cout << "G4ScoringBox::Construct() : Replicate to z direction"
             << G4endl;

    if(G4ScoringManager::GetReplicaLevel() > 2)
    {
      new G4PVReplica(elementName, fMeshElementLogical, layerLogical[1], kZAxis,
                      fNSegment[2], 2. * fSize[2] / fNSegment[2]);
    }
    else
    {
      new G4PVDivision(elementName, fMeshElementLogical, layerLogical[1],
                       kZAxis, fNSegment[2], 0.);
    }
  }
  else if(fNSegment[2] == 1)
  {
    if(verboseLevel > 9)
      G4cout << "G4ScoringBox::Construct() : Placement" << G4endl;
    new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), fMeshElementLogical,
                      elementName, layerLogical[1], false, 0);
  }
  else
    G4cerr << "ERROR : G4ScoringBox::SetupGeometry() : "
           << "invalid parameter (" << fNSegment[2] << ") "
           << "in mesh element placement." << G4endl;

  if(verboseLevel > 9)
  {
    G4cout << fSize[0] / fNSegment[0] << ", " << fSize[1] / fNSegment[1] << ", "
           << fSize[2] / fNSegment[2] << G4endl;
    G4cout << elementName << ": kZAxis, " << fNSegment[2] << ", "
           << 2. * fSize[2] / fNSegment[2] << G4endl;
  }

  // set the sensitive detector
  fMeshElementLogical->SetSensitiveDetector(fMFD);

  // vis. attributes
  auto  visatt = new G4VisAttributes(G4Colour(.5, .5, .5));
  visatt->SetVisibility(false);
  layerLogical[0]->SetVisAttributes(visatt);
  layerLogical[1]->SetVisAttributes(visatt);
  visatt->SetVisibility(true);
  fMeshElementLogical->SetVisAttributes(visatt);
}

void G4ScoringBox::List() const
{
  G4cout << "G4ScoringBox : " << fWorldName << " --- Shape: Box mesh" << G4endl;
  G4cout << " Size (x, y, z): (" << fSize[0] / cm << ", " << fSize[1] / cm
         << ", " << fSize[2] / cm << ") [cm]" << G4endl;

  G4VScoringMesh::List();
}

void G4ScoringBox::Draw(RunScore* map, G4VScoreColorMap* colorMap, G4int axflg)
{
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager != nullptr)
  {
    // cell vectors
    std::vector<std::vector<std::vector<double>>> cell;  // cell[X][Y][Z]
    std::vector<double> ez;
    for(int z = 0; z < fNSegment[2]; z++)
      ez.push_back(0.);
    std::vector<std::vector<double>> eyz;
    for(int y = 0; y < fNSegment[1]; y++)
      eyz.push_back(ez);
    for(int x = 0; x < fNSegment[0]; x++)
      cell.push_back(eyz);

    std::vector<std::vector<double>> xycell;  // xycell[X][Y]
    std::vector<double> ey;
    for(int y = 0; y < fNSegment[1]; y++)
      ey.push_back(0.);
    for(int x = 0; x < fNSegment[0]; x++)
      xycell.push_back(ey);

    std::vector<std::vector<double>> yzcell;  // yzcell[Y][Z]
    for(int y = 0; y < fNSegment[1]; y++)
      yzcell.push_back(ez);

    std::vector<std::vector<double>> xzcell;  // xzcell[X][Z]
    for(int x = 0; x < fNSegment[0]; x++)
      xzcell.push_back(ez);

    // projections
    G4int q[3];
    auto itr = map->GetMap()->begin();
    for(; itr != map->GetMap()->end(); itr++)
    {
      GetXYZ(itr->first, q);

      xycell[q[0]][q[1]] += (itr->second->sum_wx()) / fDrawUnitValue;
      yzcell[q[1]][q[2]] += (itr->second->sum_wx()) / fDrawUnitValue;
      xzcell[q[0]][q[2]] += (itr->second->sum_wx()) / fDrawUnitValue;
    }

    // search max. & min. values in each slice
    G4double xymin = DBL_MAX, yzmin = DBL_MAX, xzmin = DBL_MAX;
    G4double xymax = 0., yzmax = 0., xzmax = 0.;
    for(int x = 0; x < fNSegment[0]; x++)
    {
      for(int y = 0; y < fNSegment[1]; y++)
      {
        if(xymin > xycell[x][y])
          xymin = xycell[x][y];
        if(xymax < xycell[x][y])
          xymax = xycell[x][y];
      }
      for(int z = 0; z < fNSegment[2]; z++)
      {
        if(xzmin > xzcell[x][z])
          xzmin = xzcell[x][z];
        if(xzmax < xzcell[x][z])
          xzmax = xzcell[x][z];
      }
    }
    for(int y = 0; y < fNSegment[1]; y++)
    {
      for(int z = 0; z < fNSegment[2]; z++)
      {
        if(yzmin > yzcell[y][z])
          yzmin = yzcell[y][z];
        if(yzmax < yzcell[y][z])
          yzmax = yzcell[y][z];
      }
    }

    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);
    G4double thick = 0.01;

    G4Scale3D scale;
    if(axflg / 100 == 1)
    {
      pVisManager->BeginDraw();

      // xy plane
      if(colorMap->IfFloatMinMax())
      {
        colorMap->SetMinMax(xymin, xymax);
      }
      G4ThreeVector zhalf(0., 0., fSize[2] / fNSegment[2] - thick);
      for(int x = 0; x < fNSegment[0]; x++)
      {
        for(int y = 0; y < fNSegment[1]; y++)
        {
          G4ThreeVector pos(GetReplicaPosition(x, y, 0) - zhalf);
          G4ThreeVector pos2(GetReplicaPosition(x, y, fNSegment[2] - 1) +
                             zhalf);
          G4Transform3D trans, trans2;
          if(fRotationMatrix != nullptr)
          {
            trans = G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos);
            trans = G4Translate3D(fCenterPosition) * trans;
            trans2 =
              G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos2);
            trans2 = G4Translate3D(fCenterPosition) * trans2;
          }
          else
          {
            trans  = G4Translate3D(pos) * G4Translate3D(fCenterPosition);
            trans2 = G4Translate3D(pos2) * G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(xycell[x][y], c);
          att.SetColour(c[0], c[1], c[2]);  //, c[3]);

          G4Box xyplate("xy", fSize[0] / fNSegment[0], fSize[1] / fNSegment[1],
                        thick);
          G4Polyhedron* poly = xyplate.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(&att);
          pVisManager->Draw(*poly);

          G4Box xyplate2      = xyplate;
          G4Polyhedron* poly2 = xyplate2.GetPolyhedron();
          poly2->Transform(trans2);
          poly2->SetVisAttributes(&att);
          pVisManager->Draw(*poly2);
        }
      }
      pVisManager->EndDraw();
    }
    axflg = axflg % 100;
    if(axflg / 10 == 1)
    {
      pVisManager->BeginDraw();

      // yz plane
      if(colorMap->IfFloatMinMax())
      {
        colorMap->SetMinMax(yzmin, yzmax);
      }
      G4ThreeVector xhalf(fSize[0] / fNSegment[0] - thick, 0., 0.);
      for(int y = 0; y < fNSegment[1]; y++)
      {
        for(int z = 0; z < fNSegment[2]; z++)
        {
          G4ThreeVector pos(GetReplicaPosition(0, y, z) - xhalf);
          G4ThreeVector pos2(GetReplicaPosition(fNSegment[0] - 1, y, z) +
                             xhalf);
          G4Transform3D trans, trans2;
          if(fRotationMatrix != nullptr)
          {
            trans = G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos);
            trans = G4Translate3D(fCenterPosition) * trans;
            trans2 =
              G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos2);
            trans2 = G4Translate3D(fCenterPosition) * trans2;
          }
          else
          {
            trans  = G4Translate3D(pos) * G4Translate3D(fCenterPosition);
            trans2 = G4Translate3D(pos2) * G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(yzcell[y][z], c);
          att.SetColour(c[0], c[1], c[2]);  //, c[3]);

          G4Box yzplate("yz", thick,  // fSize[0]/fNSegment[0]*0.001,
                        fSize[1] / fNSegment[1], fSize[2] / fNSegment[2]);
          G4Polyhedron* poly = yzplate.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(&att);
          pVisManager->Draw(*poly);

          G4Box yzplate2      = yzplate;
          G4Polyhedron* poly2 = yzplate2.GetPolyhedron();
          poly2->Transform(trans2);
          poly2->SetVisAttributes(&att);
          pVisManager->Draw(*poly2);
        }
      }
      pVisManager->EndDraw();
    }
    axflg = axflg % 10;
    if(axflg == 1)
    {
      pVisManager->BeginDraw();

      // xz plane
      if(colorMap->IfFloatMinMax())
      {
        colorMap->SetMinMax(xzmin, xzmax);
      }
      G4ThreeVector yhalf(0., fSize[1] / fNSegment[1] - thick, 0.);
      for(int x = 0; x < fNSegment[0]; x++)
      {
        for(int z = 0; z < fNSegment[2]; z++)
        {
          G4ThreeVector pos(GetReplicaPosition(x, 0, z) - yhalf);
          G4ThreeVector pos2(GetReplicaPosition(x, fNSegment[1] - 1, z) +
                             yhalf);
          G4Transform3D trans, trans2;
          if(fRotationMatrix != nullptr)
          {
            trans = G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos);
            trans = G4Translate3D(fCenterPosition) * trans;
            trans2 =
              G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos2);
            trans2 = G4Translate3D(fCenterPosition) * trans2;
          }
          else
          {
            trans  = G4Translate3D(pos) * G4Translate3D(fCenterPosition);
            trans2 = G4Translate3D(pos2) * G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(xzcell[x][z], c);
          att.SetColour(c[0], c[1], c[2]);  //, c[3]);

          G4Box xzplate("xz", fSize[0] / fNSegment[0],
                        thick,  // fSize[1]/fNSegment[1]*0.001,
                        fSize[2] / fNSegment[2]);
          G4Polyhedron* poly = xzplate.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(&att);
          pVisManager->Draw(*poly);

          G4Box xzplate2      = xzplate;
          G4Polyhedron* poly2 = xzplate2.GetPolyhedron();
          poly2->Transform(trans2);
          poly2->SetVisAttributes(&att);
          pVisManager->Draw(*poly2);
        }
      }
      pVisManager->EndDraw();
    }
  }
  colorMap->SetPSUnit(fDrawUnit);
  colorMap->SetPSName(fDrawPSName);
  colorMap->DrawColorChart();
}

G4ThreeVector G4ScoringBox::GetReplicaPosition(G4int x, G4int y, G4int z)
{
  G4ThreeVector width(fSize[0] / fNSegment[0], fSize[1] / fNSegment[1],
                      fSize[2] / fNSegment[2]);
  G4ThreeVector pos(-fSize[0] + 2 * (x + 0.5) * width.x(),
                    -fSize[1] + 2 * (y + 0.5) * width.y(),
                    -fSize[2] + 2 * (z + 0.5) * width.z());

  return pos;
}

void G4ScoringBox::GetXYZ(G4int index, G4int q[3]) const
{
  q[0] = index / (fNSegment[2] * fNSegment[1]);
  q[1] = (index - q[0] * fNSegment[2] * fNSegment[1]) / fNSegment[2];
  q[2] = index - q[1] * fNSegment[2] - q[0] * fNSegment[2] * fNSegment[1];
}

G4int G4ScoringBox::GetIndex(G4int x, G4int y, G4int z) const
{
  return x + y * fNSegment[0] + z * fNSegment[0] * fNSegment[1];
}

void G4ScoringBox::DrawColumn(RunScore* map, G4VScoreColorMap* colorMap,
                              G4int idxProj, G4int idxColumn)
{
  G4int iColumn[3] = { 2, 0, 1 };
  if(idxColumn < 0 || idxColumn >= fNSegment[iColumn[idxProj]])
  {
    G4cerr << "ERROR : Column number " << idxColumn
           << " is out of scoring mesh [0," << fNSegment[iColumn[idxProj]] - 1
           << "]. Method ignored." << G4endl;
    return;
  }
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if(pVisManager != nullptr)
  {
    pVisManager->BeginDraw();

    // cell vectors
    std::vector<std::vector<std::vector<double>>> cell;  // cell[X][Y][Z]
    std::vector<double> ez;
    for(int z = 0; z < fNSegment[2]; z++)
      ez.push_back(0.);
    std::vector<std::vector<double>> eyz;
    for(int y = 0; y < fNSegment[1]; y++)
      eyz.push_back(ez);
    for(int x = 0; x < fNSegment[0]; x++)
      cell.push_back(eyz);

    std::vector<std::vector<double>> xycell;  // xycell[X][Y]
    std::vector<double> ey;
    for(int y = 0; y < fNSegment[1]; y++)
      ey.push_back(0.);
    for(int x = 0; x < fNSegment[0]; x++)
      xycell.push_back(ey);

    std::vector<std::vector<double>> yzcell;  // yzcell[Y][Z]
    for(int y = 0; y < fNSegment[1]; y++)
      yzcell.push_back(ez);

    std::vector<std::vector<double>> xzcell;  // xzcell[X][Z]
    for(int x = 0; x < fNSegment[0]; x++)
      xzcell.push_back(ez);

    // projections
    G4int q[3];
    auto itr = map->GetMap()->begin();
    for(; itr != map->GetMap()->end(); itr++)
    {
      GetXYZ(itr->first, q);

      if(idxProj == 0 && q[2] == idxColumn)
      {  // xy plane
        xycell[q[0]][q[1]] += (itr->second->sum_wx()) / fDrawUnitValue;
      }
      if(idxProj == 1 && q[0] == idxColumn)
      {  // yz plane
        yzcell[q[1]][q[2]] += (itr->second->sum_wx()) / fDrawUnitValue;
      }
      if(idxProj == 2 && q[1] == idxColumn)
      {  // zx plane
        xzcell[q[0]][q[2]] += (itr->second->sum_wx()) / fDrawUnitValue;
      }
    }

    // search max. & min. values in each slice
    G4double xymin = DBL_MAX, yzmin = DBL_MAX, xzmin = DBL_MAX;
    G4double xymax = 0., yzmax = 0., xzmax = 0.;
    for(int x = 0; x < fNSegment[0]; x++)
    {
      for(int y = 0; y < fNSegment[1]; y++)
      {
        if(xymin > xycell[x][y])
          xymin = xycell[x][y];
        if(xymax < xycell[x][y])
          xymax = xycell[x][y];
      }
      for(int z = 0; z < fNSegment[2]; z++)
      {
        if(xzmin > xzcell[x][z])
          xzmin = xzcell[x][z];
        if(xzmax < xzcell[x][z])
          xzmax = xzcell[x][z];
      }
    }
    for(int y = 0; y < fNSegment[1]; y++)
    {
      for(int z = 0; z < fNSegment[2]; z++)
      {
        if(yzmin > yzcell[y][z])
          yzmin = yzcell[y][z];
        if(yzmax < yzcell[y][z])
          yzmax = yzcell[y][z];
      }
    }

    G4VisAttributes att;
    att.SetForceSolid(true);
    att.SetForceAuxEdgeVisible(true);

    G4Scale3D scale;
    // xy plane
    if(idxProj == 0)
    {
      if(colorMap->IfFloatMinMax())
      {
        colorMap->SetMinMax(xymin, xymax);
      }
      for(int x = 0; x < fNSegment[0]; x++)
      {
        for(int y = 0; y < fNSegment[1]; y++)
        {
          G4Box xyplate("xy", fSize[0] / fNSegment[0], fSize[1] / fNSegment[1],
                        fSize[2] / fNSegment[2]);

          G4ThreeVector pos(GetReplicaPosition(x, y, idxColumn));
          G4Transform3D trans;
          if(fRotationMatrix != nullptr)
          {
            trans = G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos);
            trans = G4Translate3D(fCenterPosition) * trans;
          }
          else
          {
            trans = G4Translate3D(pos) * G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(xycell[x][y], c);
          att.SetColour(c[0], c[1], c[2]);

          G4Polyhedron* poly = xyplate.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
    }
    else
      // yz plane
      if(idxProj == 1)
    {
      if(colorMap->IfFloatMinMax())
      {
        colorMap->SetMinMax(yzmin, yzmax);
      }
      for(int y = 0; y < fNSegment[1]; y++)
      {
        for(int z = 0; z < fNSegment[2]; z++)
        {
          G4Box yzplate("yz", fSize[0] / fNSegment[0], fSize[1] / fNSegment[1],
                        fSize[2] / fNSegment[2]);

          G4ThreeVector pos(GetReplicaPosition(idxColumn, y, z));
          G4Transform3D trans;
          if(fRotationMatrix != nullptr)
          {
            trans = G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos);
            trans = G4Translate3D(fCenterPosition) * trans;
          }
          else
          {
            trans = G4Translate3D(pos) * G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(yzcell[y][z], c);
          att.SetColour(c[0], c[1], c[2]);  //, c[3]);

          G4Polyhedron* poly = yzplate.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
    }
    else
      // xz plane
      if(idxProj == 2)
    {
      if(colorMap->IfFloatMinMax())
      {
        colorMap->SetMinMax(xzmin, xzmax);
      }
      for(int x = 0; x < fNSegment[0]; x++)
      {
        for(int z = 0; z < fNSegment[2]; z++)
        {
          G4Box xzplate("xz", fSize[0] / fNSegment[0], fSize[1] / fNSegment[1],
                        fSize[2] / fNSegment[2]);

          G4ThreeVector pos(GetReplicaPosition(x, idxColumn, z));
          G4Transform3D trans;
          if(fRotationMatrix != nullptr)
          {
            trans = G4Rotate3D(*fRotationMatrix).inverse() * G4Translate3D(pos);
            trans = G4Translate3D(fCenterPosition) * trans;
          }
          else
          {
            trans = G4Translate3D(pos) * G4Translate3D(fCenterPosition);
          }
          G4double c[4];
          colorMap->GetMapColor(xzcell[x][z], c);
          att.SetColour(c[0], c[1], c[2]);  //, c[3]);

          G4Polyhedron* poly = xzplate.GetPolyhedron();
          poly->Transform(trans);
          poly->SetVisAttributes(att);
          pVisManager->Draw(*poly);
        }
      }
    }
    pVisManager->EndDraw();
  }

  colorMap->SetPSUnit(fDrawUnit);
  colorMap->SetPSName(fDrawPSName);
  colorMap->DrawColorChart();
}
