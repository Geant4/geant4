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
// $Id: PhotInConstants.hh,v 1.1 2005-11-04 13:51:36 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

