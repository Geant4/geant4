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
// $Id: RE01Trajectory.hh,v 1.1 2004/11/26 07:37:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#ifndef RE01Trajectory_h
#define RE01Trajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"     
#include "globals.hh" 
#include "G4ParticleDefinition.hh" 
#include "G4TrajectoryPoint.hh"   
#include "G4Track.hh"
#include "G4Step.hh"
#include <vector>

class G4Polyline;

typedef std::vector<G4VTrajectoryPoint*> RE01TrajectoryPointContainer;

class RE01Trajectory : public G4VTrajectory
{
 public:
   RE01Trajectory();
   RE01Trajectory(const G4Track* aTrack);
   RE01Trajectory(RE01Trajectory &);
   virtual ~RE01Trajectory();

   inline void* operator new(size_t);
   inline void  operator delete(void*);
   inline int operator == (const RE01Trajectory& right) const
   {return (this==&right);} 

   virtual G4int GetTrackID() const
   { return fTrackID; }
   virtual G4int GetParentID() const
   { return fParentID; }
   virtual G4String GetParticleName() const
   { return ParticleName; }
   virtual G4double GetCharge() const
   { return PDGCharge; }
   virtual G4int GetPDGEncoding() const
   { return PDGEncoding; }
   virtual G4ThreeVector GetInitialMomentum() const
   { return momentum; }
   virtual int GetPointEntries() const
   { return positionRecord->size(); }
   virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
   { return (*positionRecord)[i]; }

   virtual void ShowTrajectory(std::ostream& os=G4cout) const;
   virtual void DrawTrajectory(G4int i_mode=0) const;
   virtual void AppendStep(const G4Step* aStep);
   virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

   G4ParticleDefinition* GetParticleDefinition();

   inline const G4int GetTrackStatus() const
   { return fTrackStatus; }
   inline const G4ThreeVector& GetVertexPosition() const
   { return vertexPosition; }
   inline G4double GetGlobalTime() const
   { return globalTime; }
 private:
   RE01TrajectoryPointContainer* positionRecord;
   G4int                        fTrackID;
   G4int                        fParentID;
   G4int                        fTrackStatus;
   G4ParticleDefinition*        fpParticleDefinition;
   G4String                     ParticleName;
   G4double                     PDGCharge;
   G4int                        PDGEncoding;
   G4ThreeVector                momentum;
   G4ThreeVector                vertexPosition;
   G4double                     globalTime;

};

extern G4Allocator<RE01Trajectory> myTrajectoryAllocator;

inline void* RE01Trajectory::operator new(size_t)
{
  void* aTrajectory;
  aTrajectory = (void*)myTrajectoryAllocator.MallocSingle();
  return aTrajectory;
}

inline void RE01Trajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator.FreeSingle((RE01Trajectory*)aTrajectory);
}

#endif

