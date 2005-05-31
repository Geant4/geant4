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
// $Id: PhotInCalorimeterSD.cc,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#define debug
#define errors

#include "PhotInCalorimeterSD.hh"

G4int PhotInCalorimeterSD::numberOfLayers = PhotInNOfLayers;
G4int PhotInCalorimeterSD::numberOfSlabs  = PhotInNOfSlabs;

PhotInCalorimeterSD::PhotInCalorimeterSD(G4String name):G4VSensitiveDetector(name),
  SlabsCollection(0),AbsorberCollection(0),SlabsCollID(-1),AbsorberCollID(-1)
{
  for(G4int i=0; i<PhotInNumCollections; i++) collectionName.insert(PhotInCollect[i]);
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::Const: nCol="<<PhotInNumCollections<<G4endl;
#endif
}

PhotInCalorimeterSD::~PhotInCalorimeterSD()
{
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::~:S="<<SlabsCollection<<", L="<<AbsorberCollection<<G4endl;
#endif
  //delete SlabsCollection;    // Try to open later to be sure that it is not delited by G4
  //delete AbsorberCollection; // Try to open later to be sure that it is not delited by G4
}

void PhotInCalorimeterSD::Initialize(G4HCofThisEvent* HCE)
{
  // Just gives "SlabsCollection" name to the SensitiveDetectorName of the basic class
  SlabsCollection = new PhotInCalorHitsCollection(SensitiveDetectorName,collectionName[1]);
  // Why one empty hit per slab (? M.K. What for?): One Collection for all slabs!
  for(G4int i=0; i<numberOfLayers; i++)
    for(G4int j=0; j<numberOfSlabs; j++)
      SlabsCollection->insert(new PhotInCalorHit); // How to delete in the distructor?(M.K)
  if(SlabsCollID<0) SlabsCollID = GetCollectionID(0);  // If not initialized -> Initialise
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::Initialize: slabsCollectionID="<<SlabsCollID<<G4endl;
#endif
  HCE->AddHitsCollection(SlabsCollID,SlabsCollection);

  // This is not good to make the absorber active @@ in future make a flag (M.K.)
  // Just gives "AbsorberCollection" name to the SensitiveDetectorName of the basic class
  AbsorberCollection =
                    new PhotInCalorHitsCollection(SensitiveDetectorName,collectionName[0]);
  for(G4int k=0; k<numberOfLayers; k++) AbsorberCollection->insert(new PhotInCalorHit);
  if(AbsorberCollID<0) AbsorberCollID = GetCollectionID(1);// If notInitialized->Initialise
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::Initialize: layersCollectionID="<<AbsorberCollID<<G4endl;
#endif
  HCE->AddHitsCollection(AbsorberCollID, AbsorberCollection);
}

G4bool PhotInCalorimeterSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit(); // Energy deposition in the SensitiveSlab
  
  G4double stepl = aStep->GetStepLength(); // StepLength in the SensitiveSlab (can be many)
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::ProcessHits: E="<<edep<<",stepL="<<stepl<<G4endl;
#endif
  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition();
      
  // @@ G4TouchableHistoryHandler can be used
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)
                                                (aStep->GetPreStepPoint()->GetTouchable());
  //G4int LayerSlabNumber = theTouchable->GetReplicaNumber();// nL*(1+nS) (1 - absorber)
  G4int LayerNumber = theTouchable->GetReplicaNumber();    // replica # on level 0
  G4int histDepth = theTouchable->GetHistoryDepth();       // depth of enclosion
  G4int SlabNumber = -1;                                   // absorber
  if(histDepth>1) SlabNumber = theTouchable->GetReplicaNumber(1); // replica # on level 1
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::Initialize: historyDepth="<<histDepth<<", L# "<<LayerNumber
        <<", S# "<<SlabNumber<<G4endl;
#endif
		//G4int NoVolumesInLayer = 1+numberOfSlabs; // User should calculate created volumes (!)
  //G4int LayerNumber = LayerSlabNumber/NoVolumesInLayer;    // # in the AbsorberCollection
  //G4int SlabNumber  = LayerSlabNumber%NoVolumesInLayer -1; // -1 = absorber
  G4int CopyNumber = LayerNumber*numberOfSlabs+SlabNumber; // # in the SlabsCollection
  // @@ In future trac length can be defined by a flag (electrons, muons, mesons, protons)
  //G4bool charged = particle==G4Electron::ElectronDefinition()   ||
  //                 particle==G4Positron::PositronDefinition()   ||
  //                 particle==G4PionPlus::PionPlusDefinition()   ||
  //                 particle==G4PionMinus::PionMinusDefinition() ||
  //                 particle==G4MuonPlus::MuonPlusDefinition()   ||
  //                 particle==G4MuonMinus::MuonMinusDefinition() ||
  //                 particle==G4KaonPlus::KaonPlusDefinition()   ||
  //                 particle==G4KaonMinus::KaonMinusDefinition() ||
  //             				particle==G4Proton::ProtonDefinition();
  G4bool neutr   = particle==G4Neutron::NeutronDefinition();
  if(SlabNumber>-1 && LayerNumber>-1 && LayerNumber<numberOfLayers &&
     SlabNumber<numberOfSlabs) // Sensetive Slab
  {
    (*SlabsCollection)[CopyNumber]->AddEnergy(edep);
    //if(charged) (*SlabsCollection)[CopyNumber]->AddStep(stepl);
    if(neutr) (*SlabsCollection)[CopyNumber]->AddStep(stepl);
  }
  else if(SlabNumber==-1 && LayerNumber>-1 && LayerNumber<numberOfLayers) // Sens. Absorber
  {
    (*AbsorberCollection)[LayerNumber]->AddEnergy(edep);
    //if(charged) (*AbsorberCollection)[LayerNumber]->AddStep(stepl);
    if(neutr) (*AbsorberCollection)[LayerNumber]->AddStep(stepl);
  }
  else
  {
#ifdef errors
    G4cerr<<"-Warning-PhotInCalorimeterSD::ProcH: S#"<<SlabNumber<<" < "<<numberOfSlabs
          <<", L#"<<LayerNumber<<" < "<<numberOfLayers<<G4endl;
#endif
    return false;
  }
  return true;
}

void PhotInCalorimeterSD::EndOfEvent(G4HCofThisEvent*)
{
#ifdef debug
  G4cout<<"PhotInCalorimeterSD::EndOfEvent: at present it's doing nothing"<<G4endl;
#endif
  // Here something can be done in the end of the Event to summerize the information
}

void PhotInCalorimeterSD::clear()   {} // User can clean up unnecessary information 

void PhotInCalorimeterSD::DrawAll() {} // User draw the collected information 

void PhotInCalorimeterSD::PrintAll(){} // User print the collected information


