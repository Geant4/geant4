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
// $Id: PhotInLayerParameterisation.hh,v 1.3 2006/06/29 16:24:49 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
//
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal size.
//
//----------------------------------------------------------------------

#ifndef PhotInLayerParameterisation_H
#define PhotInLayerParameterisation_H 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"

#include "PhotInConstants.hh"

class PhotInLayerParameterisation : public G4VPVParameterisation
{ 
public:  // Constructors & Destructors
  PhotInLayerParameterisation();
  virtual ~PhotInLayerParameterisation(); // Means that it can be a basic class

  // -v-v-v-v-v- Virtual functions (can be overloaded) -v-v-v-v-v-v-v
  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;
  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol,
                                      const G4VTouchable* parentTouch=0);

  void SetNumberOfLayers(G4int nl)          { numberOfLayers = nl; }
  void SetAbsorberMaterial(G4Material* mat) { absMaterial = mat; }
  void SetHalfTotalThickness(G4double hz)   { hzTot = hz; }
  G4int GetNumberOfLayers()                 { return numberOfLayers; }
  G4Material* GetAbsorberMaterial()         { return absMaterial; }
  G4double GetHalfTotalThickness()          { return hzTot; }

private: // --- BODY ---
  G4int       numberOfLayers;             // Subdivision of the calorimeter in a few layers
  G4Material* absMaterial;                // Material of the absorber of the calorimeter
  G4double    hzTot;                      // TotalThickness of calorimeter to be subdivided

};


#endif


