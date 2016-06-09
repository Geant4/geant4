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
// $Id: PhotInGapParameterisation.hh,v 1.5 2006/11/22 09:44:22 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
//
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal size.
//

#ifndef PhotInGapParameterisation_H
#define PhotInGapParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

#include "PhotInConstants.hh"

class PhotInGapParameterisation : public G4VPVParameterisation
{ 
public:  // Constructors & Destructors
  PhotInGapParameterisation();
  virtual ~PhotInGapParameterisation(); // Means that it can be a basic class

  // -v-v-v-v-v- Virtual functions (can be overloaded) -v-v-v-v-v-v-v   
  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;
  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol,
                                      const G4VTouchable* parentTouch=0);

  void        SetNumberOfSlabs(G4int nl)               { numberOfSlabs = nl; }
  void        SetGapMaterial(G4Material* mat)          { gapMaterial = mat; }
  void        SetHalfTotalWidth(G4double hy)           { hyTot = hy; }
  void        SetHalfTotalLayerThickness(G4double hz)  { hzLay = hz; }
  void        SetSamplingFraction(G4double sf)         { sampFract = sf; }
  G4int       GetNumberOfSlabs() const                 { return numberOfSlabs; }
  G4Material* GetGapMaterial() const                   { return gapMaterial; }
  G4double    GetHalfTotalWidth() const                { return hyTot; }
  G4double    GetHalfTotalLayerThickness() const       { return hzLay; }
  G4double    GetSamplingFraction() const              { return sampFract; }

private: // --- BODY ---
  G4int       numberOfSlabs;           // Subdivision of the active area in a few slabs
  G4Material* gapMaterial;             // Material of the active area to be subdivided
  G4double    hyTot;                   // Total transversal halfWidth to be subdivided
  G4double    hzLay;                   // Half z-thickness of the layer of the calorimeter
  G4double    sampFract;               // Sampling fraction of the active area (in length)

};


#endif


