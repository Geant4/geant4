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
// $Id: RandomDetector.cc,v 1.2 2003-02-19 07:59:09 gcosmo Exp $
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
  G4double childRad = sqrt(3.*childDim*childDim);
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
  
  G4ThreeVector tr1( r1*cos(p1)*sin(t1), r1*sin(p1)*sin(t1), r1*cos(t1));
  G4ThreeVector tr2( r2*cos(p2)*sin(t2), r2*sin(p2)*sin(t2), r2*cos(t2));
  
  new G4PVPlacement(rm1,tr1,child1,"Child_1",aWorldLV,false,1);
  new G4PVPlacement(rm2,tr2,child2,"Child_2",aWorldLV,false,2);
  return new G4PVPlacement(0,G4ThreeVector(),aWorldLV,"Random",0,false,1);
}

