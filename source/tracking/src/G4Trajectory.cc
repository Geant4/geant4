// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Trajectory.cc,v 1.6 1999-11-11 15:38:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// ---------------------------------------------------------------
//
// G4Trajectory.cc
//
// Contact:
//   Questions and comments to this code should be sent to
//     Katsuya Amako  (e-mail: Katsuya.Amako@kek.jp)
//     Makoto  Asai   (e-mail: asai@kekvax.kek.jp)
//     Takashi Sasaki (e-mail: Takashi.Sasaki@kek.jp)
//
// ---------------------------------------------------------------


#include "G4Trajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

G4Allocator<G4Trajectory> aTrajectoryAllocator;

///////////////////////////////////////////
G4Trajectory::G4Trajectory()
///////////////////////////////////////////
{
   G4ParticleDefinition * fpParticleDefinition = 0;
   ParticleName = "";
   PDGCharge = 0;
   PDGEncoding = 0;
   fTrackID = 0;
   fParentID = 0;
   positionRecord = 0;
}

///////////////////////////////////////////
G4Trajectory::G4Trajectory(const G4Track* aTrack)
///////////////////////////////////////////
{
   G4ParticleDefinition * fpParticleDefinition = aTrack->GetDefinition();
   ParticleName = fpParticleDefinition->GetParticleName();
   PDGCharge = fpParticleDefinition->GetPDGCharge();
   PDGEncoding = fpParticleDefinition->GetPDGEncoding();
   fTrackID = aTrack->GetTrackID();
   fParentID = aTrack->GetParentID();
   positionRecord = new G4RWTPtrOrderedVector<G4VTrajectoryPoint>;
   positionRecord->insert(new G4TrajectoryPoint(aTrack->GetPosition()));
}

//////////////////////////////////////////
G4Trajectory::G4Trajectory(G4Trajectory & right)
//////////////////////////////////////////
{
  ParticleName = right.ParticleName;
  PDGCharge = right.PDGCharge;
  PDGEncoding = right.PDGEncoding;
  fTrackID = right.fTrackID;
  fParentID = right.fParentID;
  positionRecord = new G4RWTPtrOrderedVector<G4VTrajectoryPoint>;
  for(int i=0;i<right.positionRecord->entries();i++)
  {
    G4TrajectoryPoint* rightPoint = (G4TrajectoryPoint*)((*(right.positionRecord))[i]);
    positionRecord->insert(new G4TrajectoryPoint(*rightPoint));
  }
}


/////////////////////////////
G4Trajectory::~G4Trajectory()
/////////////////////////////
{
  positionRecord->clearAndDestroy();
  delete positionRecord;
}

///////////////////////////////////
void G4Trajectory::ShowTrajectory() const
///////////////////////////////////
{
   G4cout << endl << "TrackID =" << fTrackID 
        << ":ParentID=" << fParentID << endl;
   G4cout << "Particle name : " << ParticleName 
        << "  Charge : " << PDGCharge << endl;
   G4cout << "  Current trajectory has " << positionRecord->entries() 
        << " points." << endl;

   for( size_t i=0 ; i < positionRecord->entries() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
       G4cout << "Point[" << i << "]" 
            << " Position= " << aTrajectoryPoint->GetPosition() << endl;
   }
}

///////////////////////////////////////////////
void G4Trajectory::DrawTrajectory(G4int i_mode) const
///////////////////////////////////////////////
{

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   if(i_mode>=0)
   {
     G4Polyline pPolyline;
     for (int i = 0; i < positionRecord->entries() ; i++) {
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[i]);
       pos = aTrajectoryPoint->GetPosition();
       pPolyline.append( pos );
     }

     G4Colour colour;
     if(PDGCharge<0.) 
        colour = G4Colour(1.,0.,0.);
     else if(PDGCharge>0.) 
        colour = G4Colour(0.,0.,1.);
     else 
        colour = G4Colour(0.,1.,0.);

     G4VisAttributes attribs(colour);
     pPolyline.SetVisAttributes(attribs);
     if(pVVisManager) pVVisManager->Draw(pPolyline);
   }

   if(i_mode!=0)
   {
     for(int j=0; j<positionRecord->entries(); j++) {
       G4TrajectoryPoint* aTrajectoryPoint = (G4TrajectoryPoint*)((*positionRecord)[j]);
       pos = aTrajectoryPoint->GetPosition();
       G4Circle circle( pos );
       circle.SetScreenSize(0.001*i_mode);
       circle.SetFillStyle(G4Circle::filled);
       G4Colour colSpot(0.,0.,0.);
       G4VisAttributes attSpot(colSpot);
       circle.SetVisAttributes(attSpot);
       if(pVVisManager) pVVisManager->Draw(circle);
     }
   }

}

////////////////////////////////////////////
void G4Trajectory::AppendStep(const G4Step* aStep)
////////////////////////////////////////////
{
   positionRecord->append( new G4TrajectoryPoint(aStep->GetPostStepPoint()->
                                 GetPosition() ));
}
  
/////////////////////////////////////////////
G4ParticleDefinition* G4Trajectory::GetParticleDefinition()
/////////////////////////////////////////////
{
   return (G4ParticleTable::GetParticleTable()->FindParticle(ParticleName));
}

