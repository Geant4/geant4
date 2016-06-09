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
// $Id: RE01TrackInformation.hh,v 1.2 2006-06-29 17:43:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef RE01TrackInformation_h
#define RE01TrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"

class RE01TrackInformation : public G4VUserTrackInformation 
{
  public:
    RE01TrackInformation();
    RE01TrackInformation(const G4Track* aTrack);
    RE01TrackInformation(const RE01TrackInformation* aTrackInfo);
    virtual ~RE01TrackInformation();
   
    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);
    inline int operator ==(const RE01TrackInformation& right) const
    {return (this==&right);}

    RE01TrackInformation& operator =(const RE01TrackInformation& right);

    void SetSourceTrackInformation(const G4Track* aTrack);
    void Print() const;

  private:
    // Information of the primary track at the primary vertex
    G4int                 originalTrackID;  // Track ID of primary particle
    G4ParticleDefinition* particleDefinition;
    G4ThreeVector         originalPosition;
    G4ThreeVector         originalMomentum;
    G4double              originalEnergy;
    G4double              originalTime;

    G4int                 trackingStatus;
    // trackingStatus = 1 : primary or secondary track which has not yet reached to calorimeter
    //                = 0 : track which or ancester of which has reached to calorimeter

    //                = 2 : track or its ancester had once reached to calorimeter and
    //                      then escaped from it
    // Information of the track which reached to the calorimeter boundary at the boundary surface
    // This information is valid only for trackingStatus = 0 or 2
    G4int                 sourceTrackID;
    G4ParticleDefinition* sourceDefinition;
    G4ThreeVector         sourcePosition;
    G4ThreeVector         sourceMomentum;
    G4double              sourceEnergy;
    G4double              sourceTime;

  public:
    inline G4int GetOriginalTrackID() const {return originalTrackID;}
    inline G4ParticleDefinition* GetOriginalParticle() const {return particleDefinition;}
    inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
    inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
    inline G4double GetOriginalEnergy() const {return originalEnergy;}
    inline G4double GetOriginalTime() const {return originalTime;}

    inline G4int GetTrackingStatus() const {return trackingStatus;}
    inline void SetTrackingStatus(G4int i) {trackingStatus = i;}

    inline G4int GetSourceTrackID() const {return sourceTrackID;}
    inline G4ParticleDefinition* GetSourceParticle() const {return sourceDefinition;}
    inline G4ThreeVector GetSourcePosition() const {return sourcePosition;}
    inline G4ThreeVector GetSourceMomentum() const {return sourceMomentum;}
    inline G4double GetSourceEnergy() const {return sourceEnergy;}
    inline G4double GetSourceTime() const {return sourceTime;}
};

extern G4Allocator<RE01TrackInformation> aTrackInformationAllocator;

inline void* RE01TrackInformation::operator new(size_t)
{ void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void RE01TrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((RE01TrackInformation*)aTrackInfo);}

#endif

