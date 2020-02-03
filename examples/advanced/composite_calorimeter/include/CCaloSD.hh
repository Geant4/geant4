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
///////////////////////////////////////////////////////////////////////////////
// File: CCaloSD.hh
// Description: Stores hits of calorimetric type in appropriate container
//
// Use in your geometry routine: 
//    CCaloSD* caloSD = new CCaloSD(SDname, new CalorimeterOrganization);
//    G4SDManager::GetSDMpointer()->AddNewDetector(caloSD);
//  and then for every Logical Volume to be declared as sensitive
//    logVolume->SetSensitiveDetector(caloSD);
//
///////////////////////////////////////////////////////////////////////////////

#ifndef CCaloSD_h
#define CCaloSD_h 1

#include "CCalG4HitCollection.hh"
#include "CCalG4Hit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"

#include <iostream>
#include <fstream>

class CCalVOrganization;

//#define debug
 
class CCaloSD : public G4VSensitiveDetector
{

public:

  CCaloSD(G4String aSDname, CCalVOrganization* numberingScheme);
  virtual ~CCaloSD();

  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();

public:
  
  void SetPrimaryID(G4int i) {PrimaryID = i;}
  G4int  GetPrimaryID( )     {return PrimaryID;}
  void SetOrganization(CCalVOrganization* org);

private:

  G4ThreeVector SetToLocal(const G4ThreeVector& globalPoint) const;
  void getStepInfo(const G4Step* aStep);
  G4bool hitExists();
  void createNewHit();
  void updateHit();
  void StoreHit(CCalG4Hit* ahit);
  void ResetForNewPrimary();
  void summarize();
  G4double curve_LY(const G4StepPoint* stepPoint);
        
  // Data relative to primary particle (the one which triggers a shower)
  // These data are common to all Hits of a given shower.
  // One shower is made of several hits which differ by the
  // unit ID (cristal/fiber/scintillator) and the Time slice ID.

  G4ThreeVector EntrancePoint;
  G4float IncidentEnergy;
  G4int PrimID  ; //@@ ID of the primary particle.

  G4int                HCID;
  G4String             SDname;
  CCalG4HitCollection* theHC; 

  G4int              TSID; 
  CCalG4Hit*         CurrentHit;
  const G4Track*     theTrack;
  G4VPhysicalVolume* CurrentPV;
  G4VPhysicalVolume* PreviousPV;
  unsigned int       UnitID, PreviousUnitID;
  G4int              PrimaryID, TSliceID;  
  G4double           TSlice;

  const G4StepPoint*   PreStepPoint; 
  const G4StepPoint*   PostStepPoint; 
  G4float          EdepositEM, EdepositEHAD;
  G4ThreeVector  HitPoint;

  CCalVOrganization* theDescription;

};

#endif
