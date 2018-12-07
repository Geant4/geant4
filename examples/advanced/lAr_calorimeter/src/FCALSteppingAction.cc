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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <iostream>

#include "FCALSteppingAction.hh"
#include "G4SteppingManager.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "G4Event.hh"

#include "G4ThreeVector.hh"

#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALSteppingAction::FCALSteppingAction():IDold(-1),IDout(-1)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FCALSteppingAction::~FCALSteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FCALSteppingAction::UserSteppingAction(const G4Step* astep)
{ 
  // Get Edep
   G4double Edep = astep->GetTotalEnergyDeposit();

   // Get Track
     G4Track* aTrack = astep->GetTrack();

   // Get Touchable History
  G4TouchableHistory* theTouchable =  (G4TouchableHistory*)(aTrack->GetTouchable());

  // Energy deposit in FCAL1 and FCAL2
  if(Edep != 0.) 
    {
      G4VPhysicalVolume* physVol = theTouchable->GetVolume();
      
      if(strcmp(physVol->GetName(),"FCALEmModulePhysical")== 0 ||
	 strcmp(physVol->GetName(),"F1LArGapPhysical") == 0) 
	{
	  EdepFCALEm = EdepFCALEm + Edep;
	};
      
      if( (strcmp(physVol->GetName(), "FCALHadModulePhysical") == 0) ||
	  (strcmp(physVol->GetName(), "CuPlateAPhysical") == 0) ||
	  (strcmp(physVol->GetName(), "CuPlateBPhysical") == 0) ||
	  (strcmp(physVol->GetName(), "WAbsorberPhysical") == 0) ||
	  (strcmp(physVol->GetName(), "F2RodPhysical") == 0) ||
	  (strcmp(physVol->GetName(), "F2LArGapPhysical") == 0) ) 
	{
	  EdepFCALHad = EdepFCALHad + Edep;
	};
    };

  // Get Tracks properties
  G4int TrackID = aTrack->GetTrackID();
  G4int ParentID = aTrack->GetParentID();
  // Get Associated particle
  const G4DynamicParticle * aDynamicParticle = aTrack->GetDynamicParticle();
  G4ParticleDefinition * aParticle = aTrack->GetDefinition();
  G4String ParticleName = aParticle->GetParticleName();
  
  IDnow = EventNo + 10000*TrackID+ 100000000*ParentID;
  
  if(IDnow != IDold)
    {
      IDold = IDnow;
      
      // Get the primary particle
      if(TrackID==1 && ParentID==0 && (aTrack->GetCurrentStepNumber()) == 1)
	{
	  PrimaryVertex    = aTrack->GetVertexPosition(); 
	  PrimaryDirection = aTrack->GetVertexMomentumDirection();
	  
	  NSecondaries = 1;	  
	  Secondaries[NSecondaries][1] = aParticle->GetPDGEncoding();
	  Secondaries[NSecondaries][2] = PrimaryVertex.x();
	  Secondaries[NSecondaries][3] = PrimaryVertex.y();
	  Secondaries[NSecondaries][4] = PrimaryVertex.z();
	  Secondaries[NSecondaries][5] = (aDynamicParticle->GetMomentum()).x(); 
	  Secondaries[NSecondaries][6] = (aDynamicParticle->GetMomentum()).y();
	  Secondaries[NSecondaries][7] = (aDynamicParticle->GetMomentum()).z();
	  Secondaries[NSecondaries][8] = aDynamicParticle->GetTotalMomentum();
	  Secondaries[NSecondaries][9] = aDynamicParticle->GetTotalEnergy();
	  Secondaries[NSecondaries][10] = aDynamicParticle->GetKineticEnergy();
	  
	  G4cout << " ****  Primary : " << EventNo << G4endl;
	  G4cout << " Vertex : " << PrimaryVertex << G4endl;
	}
      
      
      // Get secondaries in air close to the primary tracks (DCA < 2.mm) 
      G4double DCACut = 2.*mm;
      G4String Material = aTrack->GetMaterial()->GetName();
      G4ThreeVector TrackPos = aTrack->GetVertexPosition();
      
      if(TrackID != 1 && ParentID == 1 && (strcmp(Material,"Air")==0) && (TrackPos.z() > 135.*cm)) 
	{
	  SecondaryVertex = aTrack->GetVertexPosition();
	  SecondaryDirection = aTrack->GetVertexMomentumDirection();
	  
	  // calculate DCA of secondries to primary particle
	  Distance = PrimaryVertex - SecondaryVertex ;
	  VectorProduct = PrimaryDirection.cross(SecondaryDirection);
	  if(VectorProduct == G4ThreeVector() &&  
	     PrimaryDirection != G4ThreeVector() && SecondaryDirection != G4ThreeVector()) 
	    {
	      G4ThreeVector Temp = Distance.cross(PrimaryDirection);
	      VectorProduct = Temp.cross(PrimaryDirection);
	    };	  
	  	  
VectorProductMagnitude = VectorProduct.mag();
	  if(VectorProductMagnitude == 0.) 
	    {
	      VectorProductNorm = G4ThreeVector();
	    } else {
	      VectorProductNorm = (1./VectorProduct.mag()) * VectorProduct ;
	    };	  
	  DistOfClosestApproach = Distance * VectorProductNorm ;
	  
	  if(std::abs(DistOfClosestApproach) < DCACut) 
	    {
	      NSecondaries++;	      
	      Secondaries[0][0] = NSecondaries;
	      Secondaries[NSecondaries][1] = aParticle->GetPDGEncoding();
	      Secondaries[NSecondaries][2] = (aTrack->GetVertexPosition()).x();
	      Secondaries[NSecondaries][3] = (aTrack->GetVertexPosition()).y();
	      Secondaries[NSecondaries][4] = (aTrack->GetVertexPosition()).z();
	      Secondaries[NSecondaries][5] =(aDynamicParticle->GetMomentum()).x(); 
	      Secondaries[NSecondaries][6] = (aDynamicParticle->GetMomentum()).y();
	      Secondaries[NSecondaries][7] = (aDynamicParticle->GetMomentum()).z();
	      Secondaries[NSecondaries][8] = aDynamicParticle->GetTotalMomentum();
	      Secondaries[NSecondaries][9] = aDynamicParticle->GetTotalEnergy();
	      Secondaries[NSecondaries][10] =aDynamicParticle->GetKineticEnergy();
	    };  
	};
    };


  // Get the World leaving particle
  if(aTrack->GetNextVolume() == 0) {
    if(IDnow != IDout) {
      IDout = IDnow;

      NTracks++;

      OutOfWorldTracksData[0][0] = NTracks;

      OutOfWorldTracksData[NTracks][1] = aParticle->GetPDGEncoding();

      OutOfWorldTracksData[NTracks][2] = (aTrack->GetVertexPosition()).x();
      OutOfWorldTracksData[NTracks][3] = (aTrack->GetVertexPosition()).y();
      OutOfWorldTracksData[NTracks][4] = (aTrack->GetVertexPosition()).z();

      OutOfWorldTracksData[NTracks][5] = (aDynamicParticle->GetMomentum()).x();
      OutOfWorldTracksData[NTracks][6] = (aDynamicParticle->GetMomentum()).y();
      OutOfWorldTracksData[NTracks][7] = (aDynamicParticle->GetMomentum()).z();
      
      OutOfWorldTracksData[NTracks][8] = aDynamicParticle->GetTotalMomentum();

      OutOfWorldTracksData[NTracks][9] = aDynamicParticle->GetTotalEnergy();

      OutOfWorldTracksData[NTracks][10] = aDynamicParticle->GetKineticEnergy();      
    };
  };
  
 
}

void FCALSteppingAction::initialize(G4int Nev) {
  EventNo = Nev;
  NTracks = 0;
  NSecondaries = 0;
  EdepFCALEm = EdepFCALHad = 0.;

  for(G4int i=0; i<6000; i++)
    {
      for(G4int j=0; j<11; j++) 
	{ 
	  OutOfWorldTracksData[i][j] = 0.;
	  Secondaries[i][j] = 0.; 
	}
    };
}

G4double FCALSteppingAction::GetOutOfWorldTracks(G4int i, G4int j){
  return OutOfWorldTracksData[i][j];
}

G4double FCALSteppingAction::GetSecondaries(G4int i, G4int j){
  return Secondaries[i][j];
}

G4double FCALSteppingAction::GetEdepFCAL(G4String FCAL) {
  if(strcmp(FCAL,"FCALEm") == 0) {
    return EdepFCALEm;
  } else {
    if(strcmp(FCAL,"FCALHad") == 0) {
      return EdepFCALHad;}
  }
  return 0.0; 
} 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



