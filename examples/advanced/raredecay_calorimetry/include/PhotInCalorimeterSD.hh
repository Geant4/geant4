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
// $Id: PhotInCalorimeterSD.hh,v 1.3 2006/06/29 16:24:37 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//

#ifndef PhotInCalorimeterSD_h
#define PhotInCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include "PhotInCalorHit.hh"
#include "PhotInConstants.hh"

#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"      // Headers of all Geant4 particles are here
#include "G4Track.hh"
#include "G4ios.hh"

class PhotInCalorimeterSD : public G4VSensitiveDetector
{
public:
  PhotInCalorimeterSD(G4String);
  virtual ~PhotInCalorimeterSD();

  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);
  virtual void clear();
  virtual void DrawAll();
  virtual void PrintAll();

  static void SetNumberOfLayers(G4int ln) { numberOfLayers = ln; }
  static void SetNumberOfSlabs(G4int sn)  { numberOfSlabs  = sn; }
  static G4int GetNumberOfLayers()        { return numberOfLayers; }
  static G4int GetNumberOfSlabs()         { return numberOfSlabs; }

private: // --- BODY ---
  static G4int numberOfLayers;
  static G4int numberOfSlabs;
  PhotInCalorHitsCollection*  SlabsCollection;      
  PhotInCalorHitsCollection*  AbsorberCollection;      
  G4int SlabsCollID;
  G4int AbsorberCollID;
  //G4CollectionNameVector collectionName; is from the basic class G4VSensitiveDetector
};

#endif

