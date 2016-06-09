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
// $Id: PhotInConstants.hh,v 1.2 2006/06/29 16:24:39 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//

#ifndef PhotInConstants_h
#define PhotInConstants_h 1

#include "globals.hh"

static const G4int PhotInNumCollections=2; // The # of hit collections in each section
static const G4int PhotInNumSections=3;    // The # of sections in the example
static const G4int PhotInDiNSections=PhotInNumSections*PhotInNumCollections;
static const G4int PhotInNOfLayers=10;     // Default # of Layers
static const G4int PhotInNOfSlabs=10;      // Default # of Slabs
static const G4double PhotInSampFract=0.5; // Default sampling fraction
static const G4String PhotInAbsorberName = "Absorber";
static const G4String PhotInSlabName     = "Slab";
static const G4String PhotInCalName[PhotInNumSections] = {"Sect-0","Sect-1","Sect-2"};
static const G4String PhotInRegName[PhotInNumSections] = {"Region0","Region1","Region2"};
static const G4String PhotInDetName[PhotInNumSections] = {"CalSD-0","CalSD-1","CalSD-2"};
static const G4String PhotInCollect[PhotInNumSections] =
                                                  {"AbsorberCollection","SlabsCollection"};
//@@Can make in PhotInDetectorConstruction::Construct() using PhotInDetName & PhotInCollect
static const G4String PhotInColNms[PhotInDiNSections] = {
  PhotInDetName[0]+"/"+PhotInCollect[0], PhotInDetName[0]+"/"+PhotInCollect[1],
  PhotInDetName[1]+"/"+PhotInCollect[0], PhotInDetName[1]+"/"+PhotInCollect[1],
  PhotInDetName[2]+"/"+PhotInCollect[0], PhotInDetName[2]+"/"+PhotInCollect[1]};

#endif

