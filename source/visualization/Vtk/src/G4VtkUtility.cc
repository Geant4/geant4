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
// Created by Stewart Boogert on 05/03/2023.
//

#include "G4VtkUtility.hh"

#include "G4Transform3D.hh"

#include "vtkMatrix4x4.h"
#include "vtkPlane.h"

vtkSmartPointer<vtkMatrix4x4> G4Transform3DToVtkMatrix4x4(const G4Transform3D& g4Trans)
{
  auto transform = vtkSmartPointer<vtkMatrix4x4>();

  double transformArray[16] = {g4Trans.xx(), g4Trans.xy(), g4Trans.xy(), g4Trans.dx(),
                               g4Trans.yx(), g4Trans.yy(), g4Trans.yy(), g4Trans.dy(),
                               g4Trans.zx(), g4Trans.zy(), g4Trans.zy(), g4Trans.dz(),
                               0.,           0.,           0.,           1.};
  transform->DeepCopy(transformArray);
  return transform;
}

vtkSmartPointer<vtkPlane> G4Plane3DToVtkPlane(const G4Plane3D& g4plane)
{
  auto plane = vtkSmartPointer<vtkPlane>::New();

  auto g4normal = g4plane.normal();
  auto g4point = g4plane.point();

  plane->SetNormal(g4normal.x(), g4normal.y(), g4normal.z());
  plane->SetOrigin(g4point.x(), g4point.y(), g4point.z());

  return plane;
}

G4Plane3D VtkPlaneToG4Plane3D(vtkPlane* vtkPlane)
{
  auto point = vtkPlane->GetOrigin();
  auto normal = vtkPlane->GetNormal();

  G4Plane3D plane =
    G4Plane3D(G4Normal3D(normal[0], normal[1], normal[2]), G4Point3D(point[0], point[1], point[2]));
  return plane;
}

void MaxBounds(G4double* maxBound, G4double* boundToCheck)
{
  if (boundToCheck[0] < maxBound[0]) maxBound[0] = boundToCheck[0];
  if (boundToCheck[1] > maxBound[1]) maxBound[1] = boundToCheck[1];

  if (boundToCheck[2] < maxBound[2]) maxBound[2] = boundToCheck[2];
  if (boundToCheck[3] > maxBound[3]) maxBound[3] = boundToCheck[3];

  if (boundToCheck[4] < maxBound[4]) maxBound[4] = boundToCheck[4];
  if (boundToCheck[5] > maxBound[5]) maxBound[5] = boundToCheck[5];
}
