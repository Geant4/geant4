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
// $Id: PhotInLayerParameterisation.hh,v 1.1 2005-11-04 13:51:36 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//  A parameterisation that describes a series of boxes along Z
//    The boxes have equal size.
//

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
  virtual G4Material* ComputeMaterial(const G4int copyNo, G4VPhysicalVolume* physVol);

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


