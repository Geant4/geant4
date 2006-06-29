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
// $Id: RandomDetector.cc,v 1.5 2006-06-29 17:23:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// RandomDetector
//
// Author: Martin Liendl - Martin.Liendl@cern.ch
//
// --------------------------------------------------------------
//
#include "globals.hh"
#include "Randomize.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Material.hh"

#include "RandomDetector.hh"

//RandomDetector::RandomDetector(G4int levels, G4int perLevel, G4double prop)

RandomDetector::RandomDetector(G4double prop)
 : levels_(0), 
   perLevel_(0), 
   overlapProp_(prop),
   worldDim_(10.*m)
{
}
 

RandomDetector::~RandomDetector()
{
}


G4VPhysicalVolume * RandomDetector::Construct()
{
  // Material: only one for all volumes ...
  G4double density = 1.390*g/cm3;
  G4double a = 39.95*g/mole;
  G4Material* lAr = new G4Material("liquidArgon", 18., a, density);
  
  // world volume
  G4double halfDim = worldDim_/2.;
  G4Box * aWorldBox = new G4Box("WorldBox", halfDim, halfDim, halfDim);
  G4LogicalVolume * aWorldLV = new G4LogicalVolume(aWorldBox, lAr, "WorldLV");
  
  // 2 daughters, overlapping with prop. p (protruding parent or each other)
  // G4double outer = sqrt(3.*halfDim*halfDim);
  G4double inner = halfDim;
  G4double childDim = halfDim/3.;
  G4double childRad = std::sqrt(3.*childDim*childDim);
  G4Box * aChildBox = new G4Box("ChildBox", childDim, childDim, childDim);
  G4LogicalVolume * child1 = new G4LogicalVolume(aChildBox, lAr, "Child_1_LV");
  G4LogicalVolume * child2 = new G4LogicalVolume(aChildBox, lAr, "Child_2_LV");
  
  G4bool parentOverlap = G4UniformRand() < overlapProp_ ? true : false;
  G4bool childOverlap  = G4UniformRand() < overlapProp_ ? true : false;
  
  G4ThreeVector ax1(1.,1.,1.);
  G4RotationMatrix * rm1 = new G4RotationMatrix(ax1,30.*deg);
  G4ThreeVector ax2(0.2,-1.,0.45);
  G4RotationMatrix * rm2 = new G4RotationMatrix(ax2,70.*deg);
  
  G4double t1 = G4UniformRand()*180.*deg;
  G4double p1 = G4UniformRand()*360.*deg;
  G4double t2 = G4UniformRand()*180.*deg;
  G4double p2 = G4UniformRand()*360.*deg;
  G4double r1, r2;
  
  if (parentOverlap)
  {
    r1 = inner - childRad/3.;  
  }
  else
  {
    r1 = inner - childRad - childRad/5.;
  }
  
  r2 = inner - childRad - childRad/5.;
  
  if (childOverlap)
  {
    t2 = t1;
    p2 = p1;
  }
  
  G4ThreeVector tr1( r1*std::cos(p1)*std::sin(t1), r1*std::sin(p1)*std::sin(t1), r1*std::cos(t1));
  G4ThreeVector tr2( r2*std::cos(p2)*std::sin(t2), r2*std::sin(p2)*std::sin(t2), r2*std::cos(t2));
  
  new G4PVPlacement(rm1,tr1,child1,"Child_1",aWorldLV,false,1);
  new G4PVPlacement(rm2,tr2,child2,"Child_2",aWorldLV,false,2);
  return new G4PVPlacement(0,G4ThreeVector(),aWorldLV,"Random",0,false,1);
}

