// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryTransformer.hh,v 1.1 1999-01-07 16:06:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4PromaryTransformer_h 
#define G4PromaryTransformer_h 1

#include "G4TrackVector.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"
class G4Event;
class G4PrimaryVertex;
#include "G4PrimaryParticle.hh"

class G4PrimaryTransformer
{
  public:
    G4PrimaryTransformer();
    ~G4PrimaryTransformer();
    
    G4TrackVector* GimmePrimaries(G4Event* anEvent);

  private:
    G4TrackVector TV;
    G4ParticleTable* particleTable;
    G4int verboseLevel;

  public:
    inline void SetVerboseLevel(G4int vl)
    { verboseLevel = vl; };

  private:
    void GenerateTracks(G4PrimaryVertex* primaryVertex);
    void GenerateSingleTrack(G4PrimaryParticle* primaryParticle,
              G4double x0,G4double y0,G4double z0,G4double t0);
    void SetDecayProducts(G4PrimaryParticle* mother,
                            G4DynamicParticle* motherDP);
    inline G4ParticleDefinition* GetDefinition(G4PrimaryParticle*pp)
    { 
      G4ParticleDefinition* partDef = pp->GetG4code();
      if(!partDef) partDef = particleTable->FindParticle(pp->GetPDGcode());
      return partDef;
    };
};

#endif


