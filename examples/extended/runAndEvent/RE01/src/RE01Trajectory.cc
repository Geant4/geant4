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
/// \file runAndEvent/RE01/src/RE01Trajectory.cc
/// \brief Implementation of the RE01Trajectory class
//
//
//

#include "RE01Trajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "RE01TrackInformation.hh"

G4ThreadLocal G4Allocator<RE01Trajectory> * myTrajectoryAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01Trajectory::RE01Trajectory(const G4Track* aTrack)
:G4VTrajectory(),
 fPositionRecord(0),fParticleDefinition(0)
{
   fParticleDefinition = aTrack->GetDefinition();
   fParticleName = fParticleDefinition->GetParticleName();
   fPDGCharge = fParticleDefinition->GetPDGCharge();
   fPDGEncoding = fParticleDefinition->GetPDGEncoding();
   if(fParticleName=="unknown")
   {
     G4PrimaryParticle*pp = aTrack->GetDynamicParticle()->GetPrimaryParticle();
     if(pp)
     {
       if(pp->GetCharge()<DBL_MAX) fPDGCharge = pp->GetCharge();
       fPDGEncoding = pp->GetPDGcode();
       if(pp->GetG4code()!=0)
       {
         fParticleName += " : ";
         fParticleName += pp->GetG4code()->GetParticleName();
       }
     }
   }
   fTrackID = aTrack->GetTrackID();
   RE01TrackInformation* trackInfo
    = (RE01TrackInformation*)(aTrack->GetUserInformation());
   fTrackStatus = trackInfo->GetTrackingStatus();
   if(fTrackStatus == 1)
   { fParentID = aTrack->GetParentID(); }
   else if(fTrackStatus == 2)
   { fParentID = trackInfo->GetSourceTrackID(); }
   else
   { fParentID = -1; }
   fPositionRecord = new RE01TrajectoryPointContainer();
   fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
   fMomentum = aTrack->GetMomentum();
   fVertexPosition = aTrack->GetPosition();
   fGlobalTime = aTrack->GetGlobalTime();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
RE01Trajectory::~RE01Trajectory()
{
  size_t i;
  for(i=0;i<fPositionRecord->size();i++){
    delete  (*fPositionRecord)[i];
  }
  fPositionRecord->clear();

  delete fPositionRecord;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01Trajectory::ShowTrajectory(std::ostream& os) const
{
   os << G4endl << "TrackID =" << fTrackID 
        << " : ParentID=" << fParentID << " : TrackStatus=" << fTrackStatus << G4endl;
   os << "Particle name : " << fParticleName << "  PDG code : " << fPDGEncoding
        << "  Charge : " << fPDGCharge << G4endl;
   os << "Original momentum : " <<
        G4BestUnit(fMomentum,"Energy") << G4endl;
   os << "Vertex : " << G4BestUnit(fVertexPosition,"Length")
        << "  Global time : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
   os << "  Current trajectory has " << fPositionRecord->size() 
        << " points." << G4endl;

   for( size_t i=0 ; i < fPositionRecord->size() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint = 
         (G4TrajectoryPoint*)((*fPositionRecord)[i]);
       os << "Point[" << i << "]" 
            << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01Trajectory::DrawTrajectory() const
{

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   G4Polyline pPolyline;
   for (size_t i = 0; i < fPositionRecord->size() ; i++) {
     G4TrajectoryPoint* aTrajectoryPoint = 
       (G4TrajectoryPoint*)((*fPositionRecord)[i]);
     pos = aTrajectoryPoint->GetPosition();
     pPolyline.push_back( pos );
   }

   G4Colour colour(0.2,0.2,0.2);
   if(fParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(0.,0.,1.);
   else if(fParticleDefinition==G4Electron::ElectronDefinition()
         ||fParticleDefinition==G4Positron::PositronDefinition())
      colour = G4Colour(1.,1.,0.);
   else if(fParticleDefinition==G4MuonMinus::MuonMinusDefinition()
         ||fParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      colour = G4Colour(0.,1.,0.);
   else if(fParticleDefinition->GetParticleType()=="meson")
   {
      if(fPDGCharge!=0.)
         colour = G4Colour(1.,0.,0.);
      else
         colour = G4Colour(0.5,0.,0.);
   }
   else if(fParticleDefinition->GetParticleType()=="baryon")
   {
      if(fPDGCharge!=0.)
         colour = G4Colour(0.,1.,1.);
      else
         colour = G4Colour(0.,0.5,0.5);
   }

   G4VisAttributes attribs(colour);
   pPolyline.SetVisAttributes(attribs);
   if(pVVisManager) pVVisManager->Draw(pPolyline);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
const std::map<G4String,G4AttDef>* RE01Trajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE01Trajectory",isNew);
  if (isNew) {

    G4String id("ID");
    (*store)[id] = G4AttDef(id,"Track ID","Bookkeeping","","G4int");

    G4String pid("PID");
    (*store)[pid] = G4AttDef(pid,"Parent ID","Bookkeeping","","G4int");

    G4String status("Status");
    (*store)[status] = G4AttDef(status,"Track Status","Bookkeeping","","G4int");

    G4String pn("PN");
    (*store)[pn] = G4AttDef(pn,"Particle Name","Bookkeeping","","G4String");

    G4String ch("Ch");
    (*store)[ch] = G4AttDef(ch,"Charge","Physics","e+","G4double");

    G4String pdg("PDG");
    (*store)[pdg] = G4AttDef(pdg,"PDG Encoding","Bookkeeping","","G4int");

    G4String imom("IMom");
    (*store)[imom] = G4AttDef(imom, "Momentum of track at start of trajectory",
                              "Physics","G4BestUnit","G4ThreeVector");

    G4String imag("IMag");
    (*store)[imag] = 
      G4AttDef(imag, "Magnitude of momentum of track at start of trajectory",
               "Physics","G4BestUnit","G4double");

    G4String vtxPos("VtxPos");
    (*store)[vtxPos] = G4AttDef(vtxPos, "Vertex position",
                              "Physics","G4BestUnit","G4ThreeVector");

    G4String ntp("NTP");
    (*store)[ntp] = G4AttDef(ntp,"No. of points","Bookkeeping","","G4int");

  }
  return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
std::vector<G4AttValue>* RE01Trajectory::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back
    (G4AttValue("Status",G4UIcommand::ConvertToString(fTrackStatus),""));

  values->push_back(G4AttValue("PN",fParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(fPDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(fPDGEncoding),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(fMomentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(fMomentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("VtxPos",G4BestUnit(fVertexPosition,"Length"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

  return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01Trajectory::AppendStep(const G4Step* aStep)
{
   fPositionRecord->push_back( 
       new G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition() ));
}
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void RE01Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  RE01Trajectory* seco = (RE01Trajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  //
  // initial point of the second trajectory should not be merged
  for(int i=1;i<ent;i++) 
  {
    fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
  }
  delete (*seco->fPositionRecord)[0];
  seco->fPositionRecord->clear();

}
