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

#ifndef G4VTKUTILITY_HH
#define G4VTKUTILITY_HH

#include "G4Normal3D.hh"
#include "G4Plane3D.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"
#include "G4Types.hh"

#include "vtkSmartPointer.h"

class vtkMatrix4x4;
class vtkPlane;

vtkSmartPointer<vtkMatrix4x4> G4Transform3DToVtkMatrix4x4(const G4Transform3D& g4Transformation);

vtkSmartPointer<vtkPlane> G4Plane3DToVtkPlane(const G4Plane3D& g4plane);
G4Plane3D VtkPlaneToG4Plane3D(vtkPlane* vtkPlane);
void MaxBounds(G4double* maxBound, G4double* boundToCheck);
std::size_t MakeHash(const G4ThreeVector &v);
std::size_t MakeHash(const G4Colour &c);



#endif  // G4VTKUTILITY_HH
