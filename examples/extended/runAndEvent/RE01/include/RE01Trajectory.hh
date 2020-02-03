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
/// \file runAndEvent/RE01/include/RE01Trajectory.hh
/// \brief Definition of the RE01Trajectory class
//
//
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
class G4AttDef;
class G4AttValue;

typedef std::vector<G4VTrajectoryPoint*> RE01TrajectoryPointContainer;

class RE01Trajectory : public G4VTrajectory
{
public:
  RE01Trajectory(const G4Track* aTrack);
  virtual ~RE01Trajectory();

  virtual void ShowTrajectory(std::ostream& os=G4cout) const;
  virtual void DrawTrajectory() const;
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;
  virtual void AppendStep(const G4Step* aStep);
  virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

  inline void* operator new(size_t);
  inline void  operator delete(void*);
  inline int operator == (const RE01Trajectory& right) const
  {return (this==&right);} 

  virtual G4int GetTrackID() const { return fTrackID; }
  virtual G4int GetParentID() const { return fParentID; }
  virtual G4String GetParticleName() const { return fParticleName; }
  virtual G4double GetCharge() const { return fPDGCharge; }
  virtual G4int GetPDGEncoding() const { return fPDGEncoding; }
  virtual G4ThreeVector GetInitialMomentum() const { return fMomentum; }
  virtual int GetPointEntries() const { return fPositionRecord->size(); }
  virtual G4VTrajectoryPoint* GetPoint(G4int i) const 
  { return (*fPositionRecord)[i]; }

 private:
   RE01TrajectoryPointContainer* fPositionRecord;
   G4int                        fTrackID;
   G4int                        fParentID;
   G4int                        fTrackStatus;
   G4ParticleDefinition*        fParticleDefinition;
   G4String                     fParticleName;
   G4double                     fPDGCharge;
   G4int                        fPDGEncoding;
   G4ThreeVector                fMomentum;
   G4ThreeVector                fVertexPosition;
   G4double                     fGlobalTime;

};

extern G4ThreadLocal G4Allocator<RE01Trajectory> * myTrajectoryAllocator;

inline void* RE01Trajectory::operator new(size_t)
{
  if(!myTrajectoryAllocator)
    myTrajectoryAllocator = new G4Allocator<RE01Trajectory>;
  return (void*)myTrajectoryAllocator->MallocSingle();
}

inline void RE01Trajectory::operator delete(void* aTrajectory)
{
  myTrajectoryAllocator->FreeSingle((RE01Trajectory*)aTrajectory);
}

#endif

