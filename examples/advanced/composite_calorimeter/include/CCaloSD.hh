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

#include "g4std/iostream"
#include "g4std/fstream"

class CCalVOrganization;

//#define debug
 
class CCaloSD : public G4VSensitiveDetector {

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
  
  void SetPrimaryID(int i) {PrimaryID = i;
#ifdef debug
  cout<<"CCaloSD SetPrimaryID primID ="<<i<<endl;
#endif
  }
  int  GetPrimaryID( )     {return PrimaryID;}
  void SetOrganization(CCalVOrganization* org);

private:
  
  // Data relative to primary particle (the one which triggers a shower)
  // These data are common to all Hits of a given shower.
  // One shower is made of several hits which differ by the
  // unit ID (cristal/fiber/scintillator) and the Time slice ID.

  G4ThreeVector EntrancePoint;
  float IncidentEnergy;
  G4int PrimID  ; //@@ ID of the primary particle.


private:
  G4int                HCID;
  G4String             SDname;
  CCalG4HitCollection* theHC; 

  G4int              TSID; 
  CCalG4Hit*         CurrentHit;
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
  void StoreHit(CCalG4Hit* ahit);
  void ResetForNewPrimary();
  void summarize();
      
private:
  CCalVOrganization* theDescription;

};


#endif

