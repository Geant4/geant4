///////////////////////////////////////////////////////////////////////////////
// File: G4CaloSD.hh
// Date: 11/98 
// Description: Stores hits of calorimetric type in appropriate container
//
// Use in your geometry routine: 
//    G4CaloSD* caloSD = new G4CaloSD(SDname, new CalorimeterOrganization);
//    G4SDManager::GetSDMpointer()->AddNewDetector(caloSD);
//  and then for every Logical Volume to be declared as sensitive
//    logVolume->SetSensitiveDetector(caloSD);
//
// Modifications: 27/04/00 SB 
///////////////////////////////////////////////////////////////////////////////

#ifndef G4CaloSD_h
#define G4CaloSD_h 1

#include "G4CaloHitsCollection.hh"
#include "G4CaloHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"

#include <iostream>
#include <fstream>

class VDetectorOrganization;

//#define debug
 
class G4CaloSD : public G4VSensitiveDetector {

public:
  G4CaloSD(G4String aSDname, VDetectorOrganization* numberingScheme);
  virtual ~G4CaloSD();

  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();

public:
  
  void SetPrimaryID(int i) {PrimaryID = i;
#ifdef debug
  cout<<"G4CaloSD SetPrimaryID primID ="<<i<<endl;
#endif
  }
  int  GetPrimaryID( )     {return PrimaryID;}
  void SetOrganization(VDetectorOrganization* org);

private:
  
  // Data relative to primary particle (the one which triggers a shower)
  // These data are common to all Hits of a given shower.
  // One shower is made of several hits which differ by the
  // unit ID (cristal/fiber/scintillator) and the Time slice ID.

  G4ThreeVector EntrancePoint;
  float IncidentEnergy;
  G4int PrimID  ; //@@ ID of the primary particle.


private:
  G4int                 HCID;
  G4String              SDname;
  G4CaloHitsCollection* theHC; 

  G4int              TSID; 
  G4CaloHit*         CurrentHit;
  G4Track*           theTrack;
  G4VPhysicalVolume* CurrentPV;
  G4VPhysicalVolume* PreviousPV;
  unsigned int       UnitID, PreviousUnitID;
  G4int              PrimaryID, TSliceID;  
  G4double           TSlice;

  G4StepPoint*   PreStepPoint  ; 
  G4StepPoint*   PostStepPoint ; 
  float          EdepositEM, EdepositEHAD     ;
  G4ThreeVector  HitPoint      ;	 
 
private:

  G4ThreeVector SetToLocal(G4ThreeVector globalPoint);
  void getStepInfo(G4Step* aStep);
  G4bool hitExists();
  void createNewHit();
  void updateHit();
  void StoreHit(G4CaloHit* ahit);
  void ResetForNewPrimary();
  void summarize();
      
private:
  VDetectorOrganization* theDescription;

};


#endif

