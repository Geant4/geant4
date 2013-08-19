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
// $Id: G4MagneticFieldModel.cc 66839 2013-01-14 13:19:31Z allison $
//
// 
// John Allison  17th August 2013
// Model that knows how to draw the magnetic field.

#include "G4MagneticFieldModel.hh"

#include "G4PhysicalConstants.hh"
#include "G4VGraphicsScene.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4Polyhedron.hh"
#include "G4Vector3D.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"

G4MagneticFieldModel::~G4MagneticFieldModel ()
{
}

G4MagneticFieldModel::G4MagneticFieldModel ()
{
  fType = "G4MagneticFieldModel";
  fGlobalTag = fType;
  fGlobalDescription = fType;
}

void G4MagneticFieldModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler)
{
  const G4VisExtent& extent = sceneHandler.GetExtent();
  G4cout
  << "G4MagneticFieldModel::DescribeYourselfTo: Scene extent: "
  << extent
  << G4endl;
  sceneHandler.BeginPrimitives();
  sceneHandler.EndPrimitives();
}
