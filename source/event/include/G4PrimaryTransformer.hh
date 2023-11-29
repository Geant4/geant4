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
// G4PrimaryTransformer
//
// Class description:
//
// This class is exclusively used by G4EventManager for the conversion
// from G4PrimaryVertex/G4PrimaryParticle to G4DynamicParticle/G4Track.

// Author: Makoto Asai, 1999
// --------------------------------------------------------------------
#ifndef G4PrimaryTransformer_hh 
#define G4PrimaryTransformer_hh 1

#include "globals.hh"
#include "G4TrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4DynamicParticle.hh"

class G4Event;
class G4PrimaryVertex;

class G4PrimaryTransformer
{
  public:

    G4PrimaryTransformer();
    virtual ~G4PrimaryTransformer() = default;
    
    G4TrackVector* GimmePrimaries(G4Event* anEvent, G4int trackIDCounter=0);
    void CheckUnknown();

    inline void SetVerboseLevel(G4int vl)
      { verboseLevel = vl; }

    void SetUnknnownParticleDefined(G4bool vl);
      // By invoking this method, the user can alter the treatment of unknown
      // particles. The ideal place to invoke this method is in the
      // BeginOfRunAction().

    inline G4bool GetUnknownParticleDefined() const
      { return unknownParticleDefined; }

  protected:

    void GenerateTracks(G4PrimaryVertex* primaryVertex);
    void GenerateSingleTrack(G4PrimaryParticle* primaryParticle,
                             G4double x0, G4double y0, G4double z0,
                             G4double t0, G4double wv);
    void SetDecayProducts(G4PrimaryParticle* mother,
                          G4DynamicParticle* motherDP);
    G4bool CheckDynamicParticle(G4DynamicParticle*DP);

    // Following two virtual methods are provided to customize the use
    // of G4PrimaryTransformer for particle types exotic to Geant4.

    virtual G4ParticleDefinition* GetDefinition(G4PrimaryParticle* pp);
      // Return appropriate G4ParticleDefinition w.r.t. the primary particle. 
      // If nullptr is returned, the particle will not be treated as a track,
      // but its daughters will be examined in case it has "pre-assigned
      // decay products".

    virtual G4bool IsGoodForTrack(G4ParticleDefinition* pd);
      // Return true if a primary particle should be converted into a track.
      // By default, all particles of non-shortlived and shortlived with
      // valid decay tables are converted.

  protected:

    G4TrackVector TV;
    G4ParticleTable* particleTable = nullptr;

    G4ParticleDefinition* unknown = nullptr;
    G4ParticleDefinition* opticalphoton = nullptr;
    G4int verboseLevel = 0;
    G4int trackID = 0;
    G4int nWarn = 0;
    G4bool unknownParticleDefined = false;
    G4bool opticalphotonDefined = false;
};

#endif
