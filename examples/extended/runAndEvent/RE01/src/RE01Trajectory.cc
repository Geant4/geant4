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
// $Id: RE01Trajectory.cc,v 1.5 2010-11-08 22:16:28 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

G4Allocator<RE01Trajectory> myTrajectoryAllocator;

RE01Trajectory::RE01Trajectory()
:G4VTrajectory()
{
   fpParticleDefinition = 0;
   ParticleName = "";
   PDGCharge = 0;
   PDGEncoding = 0;
   fTrackID = 0;
   fParentID = 0;
   fTrackStatus = 0;
   positionRecord = 0;
   momentum = G4ThreeVector(0.,0.,0.);
   vertexPosition = G4ThreeVector(0.,0.,0.);
   globalTime = 0.;
}

RE01Trajectory::RE01Trajectory(const G4Track* aTrack)
:G4VTrajectory()
{
   fpParticleDefinition = aTrack->GetDefinition();
   ParticleName = fpParticleDefinition->GetParticleName();
   PDGCharge = fpParticleDefinition->GetPDGCharge();
   PDGEncoding = fpParticleDefinition->GetPDGEncoding();
   if(ParticleName=="unknown")
   {
     G4PrimaryParticle*pp = aTrack->GetDynamicParticle()->GetPrimaryParticle();
     if(pp)
     {
       if(pp->GetCharge()<DBL_MAX) PDGCharge = pp->GetCharge();
       PDGEncoding = pp->GetPDGcode();
       if(pp->GetG4code()!=0)
       {
         ParticleName += " : ";
         ParticleName += pp->GetG4code()->GetParticleName();
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
   positionRecord = new RE01TrajectoryPointContainer();
   positionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
   momentum = aTrack->GetMomentum();
   vertexPosition = aTrack->GetPosition();
   globalTime = aTrack->GetGlobalTime();
}

RE01Trajectory::RE01Trajectory(RE01Trajectory & right)
:G4VTrajectory()
{
  ParticleName = right.ParticleName;
  fpParticleDefinition = right.fpParticleDefinition;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  fTrackStatus = right.fTrackStatus;
  positionRecord = new RE01TrajectoryPointContainer();
  for(size_t i=0;i<right.positionRecord->size();i++)
  {
    G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.positionRecord))[i]);
    positionRecord->push_back(new G4TrajectoryPoint(*rightPoint));
  }
   momentum = right.momentum;
   vertexPosition = right.vertexPosition;
   globalTime = right.globalTime;
}

RE01Trajectory::~RE01Trajectory()
{
  size_t i;
  for(i=0;i<positionRecord->size();i++){
    delete  (*positionRecord)[i];
  }
  positionRecord->clear();

  delete positionRecord;
}

void RE01Trajectory::ShowTrajectory(std::ostream& os) const
{
   os << G4endl << "TrackID =" << fTrackID 
        << " : ParentID=" << fParentID << " : TrackStatus=" << fTrackStatus << G4endl;
   os << "Particle name : " << ParticleName << "  PDG code : " << PDGEncoding
        << "  Charge : " << PDGCharge << G4endl;
   os << "Original momentum : " <<
        G4BestUnit(momentum,"Energy") << G4endl;
   os << "Vertex : " << G4BestUnit(vertexPosition,"Length")
        << "  Global time : " << G4BestUnit(globalTime,"Time") << G4endl;
   os << "  Current trajectory has " << positionRecord->size() 
        << " points." << G4endl;

   for( size_t i=0 ; i < positionRecord->size() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
       os << "Point[" << i << "]" 
            << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   }
}

void RE01Trajectory::DrawTrajectory(G4int /*i_mode*/) const
{

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   G4Polyline pPolyline;
   for (size_t i = 0; i < positionRecord->size() ; i++) {
     G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
     pos = aTrajectoryPoint->GetPosition();
     pPolyline.push_back( pos );
   }

   G4Colour colour(0.2,0.2,0.2);
   if(fpParticleDefinition==G4Gamma::GammaDefinition())
      colour = G4Colour(0.,0.,1.);
   else if(fpParticleDefinition==G4Electron::ElectronDefinition()
         ||fpParticleDefinition==G4Positron::PositronDefinition())
      colour = G4Colour(1.,1.,0.);
   else if(fpParticleDefinition==G4MuonMinus::MuonMinusDefinition()
         ||fpParticleDefinition==G4MuonPlus::MuonPlusDefinition())
      colour = G4Colour(0.,1.,0.);
   else if(fpParticleDefinition->GetParticleType()=="meson")
   {
      if(PDGCharge!=0.)
         colour = G4Colour(1.,0.,0.);
      else
         colour = G4Colour(0.5,0.,0.);
   }
   else if(fpParticleDefinition->GetParticleType()=="baryon")
   {
      if(PDGCharge!=0.)
         colour = G4Colour(0.,1.,1.);
      else
         colour = G4Colour(0.,0.5,0.5);
   }

   G4VisAttributes attribs(colour);
   pPolyline.SetVisAttributes(attribs);
   if(pVVisManager) pVVisManager->Draw(pPolyline);
}

const std::map<G4String,G4AttDef>* RE01Trajectory::GetAttDefs() const
{
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("RE01Trajectory",isNew);
  if (isNew) {

    G4String ID("ID");
    (*store)[ID] = G4AttDef(ID,"Track ID","Bookkeeping","","G4int");

    G4String PID("PID");
    (*store)[PID] = G4AttDef(PID,"Parent ID","Bookkeeping","","G4int");

    G4String Status("Status");
    (*store)[Status] = G4AttDef(Status,"Track Status","Bookkeeping","","G4int");

    G4String PN("PN");
    (*store)[PN] = G4AttDef(PN,"Particle Name","Bookkeeping","","G4String");

    G4String Ch("Ch");
    (*store)[Ch] = G4AttDef(Ch,"Charge","Physics","e+","G4double");

    G4String PDG("PDG");
    (*store)[PDG] = G4AttDef(PDG,"PDG Encoding","Bookkeeping","","G4int");

    G4String IMom("IMom");
    (*store)[IMom] = G4AttDef(IMom, "Momentum of track at start of trajectory",
			      "Physics","G4BestUnit","G4ThreeVector");

    G4String IMag("IMag");
    (*store)[IMag] = 
      G4AttDef(IMag, "Magnitude of momentum of track at start of trajectory",
	       "Physics","G4BestUnit","G4double");

    G4String VtxPos("VtxPos");
    (*store)[VtxPos] = G4AttDef(VtxPos, "Vertex position",
			      "Physics","G4BestUnit","G4ThreeVector");

    G4String NTP("NTP");
    (*store)[NTP] = G4AttDef(NTP,"No. of points","Bookkeeping","","G4int");

  }
  return store;
}

std::vector<G4AttValue>* RE01Trajectory::CreateAttValues() const
{
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back
    (G4AttValue("Status",G4UIcommand::ConvertToString(fTrackStatus),""));

  values->push_back(G4AttValue("PN",ParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(PDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(PDGEncoding),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(momentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(momentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("VtxPos",G4BestUnit(vertexPosition,"Length"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

  return values;
}

void RE01Trajectory::AppendStep(const G4Step* aStep)
{
   positionRecord->push_back( new G4TrajectoryPoint(aStep->GetPostStepPoint()->
                                 GetPosition() ));
}
  
G4ParticleDefinition* RE01Trajectory::GetParticleDefinition()
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

void RE01Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory)
{
  if(!secondTrajectory) return;

  RE01Trajectory* seco = (RE01Trajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  for(int i=1;i<ent;i++) // initial point of the second trajectory should not be merged
  {
    positionRecord->push_back((*(seco->positionRecord))[i]);
  }
  delete (*seco->positionRecord)[0];
  seco->positionRecord->clear();

}



