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
// $Id: HadrontherapyCalorimeterSD.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#include "HadrontherapyCalorimeterSD.hh"
#include "G4HCofThisEvent.hh"
#include "HadrontherapyEventAction.hh"
#include "HadrontherapyHit.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "HadrontherapyHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4Positron.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "HadrontherapyRunAction.hh"  
#include "G4ios.hh"

HadrontherapyCalorimeterSD::HadrontherapyCalorimeterSD(G4String name, 
						       HadrontherapyDetectorConstruction* det)
:G4VSensitiveDetector(name),Detector(det),
delta(1.0e-6*mm),evno(0),evnOld(-1),trIDold(-1),NbOfLayer(20000) 
{
  collectionName.insert("CalCollection");
  
  p_Run = new HadrontherapyRunAction();
}

// -------------------------------------------------------------------
HadrontherapyCalorimeterSD::~HadrontherapyCalorimeterSD()
{
  delete  p_Run;
}

// -----------------------------------------------------------------
void HadrontherapyCalorimeterSD::Initialize(G4HCofThisEvent*)
{
  CalCollection = new HadrontherapyHitsCollection
    (SensitiveDetectorName,collectionName[0]); 
  
  for (G4int i=0; i < NbOfLayer; i++) 
    {
      energy[i]  = 0.0;
    }
  
for (G4int k=0; k < NbOfLayer; k++) 
    {
      sliceID[k] = -1;
    }
  
depthMax = (0.002*NbOfLayer - delta); // 0.2 is the actual dosimeter thickness
}

// ---------------------------------------------------------------------
G4bool HadrontherapyCalorimeterSD::ProcessHits(G4Step* aStep,G4TouchableHistory*) 
{
  
G4int j = aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetCopyNo();  
G4double edep = aStep->GetTotalEnergyDeposit();
 if (edep ==0) return false;
 if(sliceID[j]== -1)
   {
     HadrontherapyHit* newHit = new HadrontherapyHit(j);
     
     newHit->SetEdep(edep);
     G4int barbatrucco_ID = CalCollection->insert(newHit);
     sliceID[j] = barbatrucco_ID - 1;
     
   }
 else
   {
     (*CalCollection)[sliceID[j]]->AddEdep(edep);
   }
   
 if (j < 0 || j >= NbOfLayer)
   
   {
     G4cout << " Warning!! HadrontherapyCalorimeterSD cannot add " 
	    << edep 
	    << " to the layer Nb " << j << G4endl;
   }
 else
   
   {
     energy[j] = energy[j] +  edep;
   }
 
const G4Track* track = aStep->GetTrack();

G4int trIDnow  = track->GetTrackID();
G4double tkin  = track->GetKineticEnergy(); 
G4double xend  = aStep->GetPostStepPoint()->GetPosition().x();

G4bool stop = false;
G4bool primary = false;
G4bool outAbs = false;

  if(track->GetNextVolume()->GetName() == "TreatmentRoom"
     && (xend <= delta || xend >= depthMax-delta)  ) outAbs = true;

  if(tkin == 0.0) stop = true;

 // Here checks if the particle is primary

  if(aStep->GetTrack()->GetParentID() == 0)
    {
      primary = true;
    }
  
 // new particle
  
  if(trIDnow != trIDold || evno != evnOld) 
    {
      trIDold = trIDnow;
      evnOld  = evno;
      part_is_out = true;
      tkinold = aStep->GetPreStepPoint()->GetKineticEnergy();
    }
  
 // Primary particle stop 
  if(primary && (stop || outAbs)) {
  }

  if(outAbs && part_is_out) {
    G4double e = tkin;
    if(track->GetDefinition() == G4Positron::Positron()) 
      e += 2.*electron_mass_c2;
    if(xend > depthMax-delta)       leakEnergy += e;
    else if(xend < delta)       backEnergy += e;
    part_is_out = false;
  }

  return true;
}

// ----------------------------------------------------------------------
void HadrontherapyCalorimeterSD::EndOfEvent(G4HCofThisEvent* HCE)
{

G4int i;
G4double etot = 0.0;
 if (NbOfLayer > 0)
   {
     G4double s = 0.002;  //Actual thickness in mm of the dosimeter
     G4double z = -0.5 * s;
     
     for (i = 0; i < NbOfLayer; i ++)
       {
	 z += s;
	 etot += energy[i];
       }
   }
  
  //*******************************
  // Total energy of a particle
  //  that reach the dosemeter
  //*******************************
 
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = GetCollectionID(0); 
    }
  HCE->AddHitsCollection(HCID, CalCollection);
  
} 

// -------------------------------------------------------------------------

void HadrontherapyCalorimeterSD::clear()
{} 

void HadrontherapyCalorimeterSD::PrintAll()
{} 




