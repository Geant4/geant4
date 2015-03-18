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

#include "Test2PhantomSD.hh"
#include "Test2PhantomHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PSDirectionFlag.hh"


Test2PhantomSD::Test2PhantomSD(G4String name, G4int segment[3])
  :G4VSensitiveDetector(name) {

  nSegment[0]=segment[0];
  nSegment[1]=segment[1];
  nSegment[2]=segment[2];
  G4String HCname;
  collectionName.insert(HCname = "PhantomCollection");

}

Test2PhantomSD::~Test2PhantomSD() {
  ;
}

void Test2PhantomSD::Initialize(G4HCofThisEvent *) {

  fPhantomCollection = new Test2PhantomHitsCollection(SensitiveDetectorName,
						      collectionName[0]); 
  verboseLevel = 0;

}

G4bool Test2PhantomSD::ProcessHits(G4Step * aStep, G4TouchableHistory *) {

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double trklen = aStep->GetStepLength();
  if(edep > 0. || trklen > 0.) {
    edep /= MeV;
    trklen /= mm;
    // total energy deposit
    G4double weight = aStep->GetPreStepPoint()->GetWeight();
    G4double density = aStep->GetPreStepPoint()->GetMaterial()->GetDensity();
    G4double volume = GetVolume(aStep);
    G4double area   = GetArea(aStep);
    G4double charge = aStep->GetPreStepPoint()->GetCharge();
    
    G4double dose = edep/(density*volume)/gray;

    G4double flatSurfCurr = 0.0;    
    G4double flatSurfFlux= 0.0;
    G4int dirFlag = IsSelectedSurface(aStep);
    if ( dirFlag > 0 ){
      flatSurfCurr = 1.0*weight/area/(1./cm2);
      G4double anglefac = GetAngleFactor(aStep,dirFlag);
      flatSurfFlux = weight/anglefac/area/(1./cm2);
    }

    G4double passageCellCurr = 0.0;
    G4double passageCellFlux = 0.0;
    if ( IsPassed(aStep) ){
      passageCellCurr = weight;
      passageCellFlux = fCellTrack/volume/(1./cm2); // fCellTrack is calculated in IsPassed().
    }

    G4double cellFlux = trklen*weight/volume/(1./cm2);

    G4double nOfSecondary=0.0;
    if ( IsSecondary(aStep) ){
      nOfSecondary = weight;
    }

    G4double cellCharge = 0.0;
    if ( IsEnterOrFirstStep(aStep) ){
      cellCharge = charge*weight;
    }else if ( IsExit(aStep) ){
      cellCharge = -1.*charge*weight;
    }

      if(verboseLevel > 1) G4cout << "Next step edep [MeV] = " << edep/MeV << G4endl;

    G4TouchableHistory * hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
    G4int copyIDinX = hist->GetReplicaNumber(2);
    G4int copyIDinY = hist->GetReplicaNumber(1);
    G4int copyIDinZ = hist->GetReplicaNumber(0);

    Test2PhantomHit* phantomHit
      = new Test2PhantomHit(copyIDinX, copyIDinY, copyIDinZ,nSegment);
    phantomHit->SetEdep(edep);
    phantomHit->SetTrackLength(trklen);
    phantomHit->SetParticleName(aStep->GetTrack()->GetParticleDefinition()->GetParticleName());
    phantomHit->SetDose(dose);

    phantomHit->SetFlatSurfaceCurrent(flatSurfCurr);
    phantomHit->SetPassageCellCurrent(passageCellCurr);
    phantomHit->SetFlatSurfaceFlux(flatSurfFlux);
    phantomHit->SetCellFlux(cellFlux);
    phantomHit->SetPassageCellFlux(passageCellFlux);
    phantomHit->SetNofSecondary(nOfSecondary);
    phantomHit->SetCellCharge(cellCharge);

    fPhantomCollection->insert(phantomHit);
    if(verboseLevel > 0) {
      G4cout << " A hit of Test2PhantomHit is created in copy-id (" 
	     << copyIDinX << ", " << copyIDinY << ", " << copyIDinZ
	     << ") at " << GetFullPathName() << G4endl;
    }
  }

  return true;
}

void Test2PhantomSD::EndOfEvent(G4HCofThisEvent * HCE) {

  static G4int HCID = -1;
  if(HCID < 0) { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection( HCID, fPhantomCollection );
}

void Test2PhantomSD::clear() {
  ;
} 

void Test2PhantomSD::DrawAll() {
  ;
} 

void Test2PhantomSD::PrintAll() {
  ;
} 


G4VSolid* Test2PhantomSD::GetSolid(G4Step* aStep){
  G4VPhysicalVolume* physVol = aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid = 0;
  G4int indexDepth = 0;
  if(physParam)
    { // for parameterized volume                                                 
      G4int idx = ((G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable()))
	->GetReplicaNumber(indexDepth);
      if(idx<0)
	{
	  G4Exception("G4PSDoseDeposit","G4PSDoseDeposit::ProcessHits",JustWarning,
		      "Incorrect replica number");
	  G4cerr << " --- GetReplicaNumber : " << idx << G4endl;
	}
      solid = physParam->ComputeSolid(idx, physVol);
      solid->ComputeDimensions(physParam,idx,physVol);
    }
  else
    { // for ordinary volume                                                      
      solid = physVol->GetLogicalVolume()->GetSolid();
    }

  return solid;
}

G4double Test2PhantomSD::GetVolume(G4Step* aStep){
  G4VSolid* solid = GetSolid(aStep);
  return (solid->GetCubicVolume());
}

G4double Test2PhantomSD::GetArea(G4Step* aStep){
  G4Box* boxSolid = (G4Box*)(GetSolid(aStep));
  return (4.*boxSolid->GetXHalfLength()*boxSolid->GetYHalfLength());
}

G4int Test2PhantomSD::IsSelectedSurface(G4Step* aStep){
  G4Box* boxSolid = (G4Box*)GetSolid(aStep);
  G4TouchableHandle theTouchable =
    aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Entering Geometry                                                        
    G4ThreeVector stppos1= aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 =
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    if(std::fabs( localpos1.z() + boxSolid->GetZHalfLength())<kCarTolerance ){
      return fCurrent_In;
    }
  }

  if (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary ){
    // Exiting Geometry                                                         
    G4ThreeVector stppos2= aStep->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos2 =
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos2);
    if(std::fabs( localpos2.z() + boxSolid->GetZHalfLength())<kCarTolerance ){
      return fCurrent_Out;
    }
  }

  return -1;
}

G4bool Test2PhantomSD::IsPassed(G4Step* aStep){
  G4bool Passed = FALSE;

  G4bool IsEnter = aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary;
  G4bool IsExit  = aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary;
  G4double trklength  = aStep->GetStepLength();
  trklength *= aStep->GetPreStepPoint()->GetWeight();
  G4int  trkid  = aStep->GetTrack()->GetTrackID();

  if ( IsEnter &&IsExit ){         // Passed at one step                        
    fCellTrack = trklength;
    Passed = TRUE;
  }else if ( IsEnter ){            // Enter a new geometry                      
    fCurrentTrkID = trkid;         // Resetting the current track.              
    fCellTrack = trklength;
  }else if ( IsExit ){             // Exit a current geometry                   
    if ( fCurrentTrkID == trkid ) {
      fCellTrack += trklength;
      Passed = TRUE;               // if the track is same as entered.          
    }
  }else{                           // Inside geometry                           
    if ( fCurrentTrkID == trkid ){ // Adding the track length to current one ,  
      fCellTrack = trklength;
    }
  }
  return Passed;
}


G4double Test2PhantomSD::GetAngleFactor(G4Step* aStep, G4int dirFlag){
  G4StepPoint* thisStep=0;
  if ( dirFlag == fFlux_In ){
    thisStep = aStep->GetPreStepPoint();
  }else if ( dirFlag == fFlux_Out ){
    thisStep = aStep->GetPostStepPoint();
  }else{
    return FALSE;
  }
  G4TouchableHandle theTouchable = thisStep->GetTouchableHandle();
  G4ThreeVector pdirection = thisStep->GetMomentumDirection();
  G4ThreeVector localdir  =
    theTouchable->GetHistory()->GetTopTransform().TransformAxis(pdirection);
  //                                                                        
  G4double angleFactor = localdir.z();
  if ( angleFactor < 0 ) angleFactor *= -1.;
  return angleFactor;
}

G4bool Test2PhantomSD::IsSecondary(G4Step* aStep){
  if ( aStep->GetTrack()->GetCurrentStepNumber() != 1) return FALSE;
  //- check for this is not a primary particle. e.g. ParentID > 0 .             
  if ( aStep->GetTrack()->GetParentID() == 0 ) return FALSE;
  //- check the particle if the partifle definition is given.                   

  return TRUE;
}
G4bool Test2PhantomSD::IsEnterOrFirstStep(G4Step* aStep){
  // Enter or First step of primary.                                          
  return (aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary
      || ( aStep->GetTrack()->GetParentID() == 0 &&
	   aStep->GetTrack()->GetCurrentStepNumber() == 1 ) );
}  
G4bool Test2PhantomSD::IsExit(G4Step* aStep){
  return (aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
}
