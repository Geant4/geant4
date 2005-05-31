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
// $Id: PhotInGapParameterisation.hh,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol);

  void        SetNumberOfSlabs(G4int nl)               { numberOfSlabs = nl; }
  void        SetGapMaterial(G4Material* mat)          { gapMaterial = mat; }
  void        SetHalfTotalWidth(G4double hy)           { hyTot = hy; }
  void        SetHalfTotalLayerThickness(G4double hz)  { hzLay = hz; }
  void        SetSamplingFraction(G4double sf)         { sampFract = sf; }
  G4int       GetNumberOfSlabs()                       { return numberOfSlabs; }
  G4Material* GetGapMaterial()                         { return gapMaterial; }
  G4double    GetHalfTotalWidth(G4double hy)           { return hyTot; }
  G4double    GetHalfTotalLayerThickness(G4double hz)  { return hzLay; }
  G4double    GetSamplingFraction(G4double sf)         { return sampFract; }

private: // --- BODY ---
  G4int       numberOfSlabs;           // Subdivision of the active area in a few slabs
  G4Material* gapMaterial;             // Material of the active area to be subdivided
  G4double    hyTot;                   // Total transversal halfWidth to be subdivided
  G4double    hzLay;                   // Half z-thickness of the layer of the calorimeter
  G4double    sampFract;               // Sampling fraction of the active area (in length)

};


#endif


