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
// $Id: PhotInCalorimeterSD.hh,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

