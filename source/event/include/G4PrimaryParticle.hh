// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryParticle.hh,v 1.1 1999-01-07 16:06:33 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//


#ifndef G4PrimaryParticle_h
#define G4PrimaryParticle_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4ParticleDefinition;

class G4PrimaryParticle 
{
  public:
      inline void *operator new(size_t);
      inline void operator delete(void *aStackedTrack);

      G4PrimaryParticle();
      G4PrimaryParticle(G4int Pcode);
      G4PrimaryParticle(G4int Pcode,
                        G4double px,G4double py,G4double pz);
      G4PrimaryParticle(G4ParticleDefinition* Gcode);
      G4PrimaryParticle(G4ParticleDefinition* Gcode,
                        G4double px,G4double py,G4double pz);
      ~G4PrimaryParticle();

      const G4PrimaryParticle & operator=(const G4PrimaryParticle &right);
      int operator==(const G4PrimaryParticle &right) const;
      int operator!=(const G4PrimaryParticle &right) const;

      void Print() const;

  private:
      G4int PDGcode;
      G4ParticleDefinition * G4code;
      G4double Px;
      G4double Py;
      G4double Pz;
      
      G4PrimaryParticle * nextParticle;
      G4PrimaryParticle * daughterParticle;

      G4int trackID;  // This will be set if this particle is
                      // sent to G4EventManager and converted to
                      // G4Track. Otherwise = -1.
      G4double mass;  // This is just for book keeping.
                      // This will not be used but the mass in
                      // G4ParticleDefinition will be used.

      G4double polX;
      G4double polY;
      G4double polZ;

  public:
      void SetPDGcode(G4int Pcode);
      inline G4int GetPDGcode() const
      { return PDGcode; };
      void SetG4code(G4ParticleDefinition * Gcode);
      inline G4ParticleDefinition * GetG4code() const
      { return G4code; };
      inline void SetMomentum(G4double px, G4double py, G4double pz)
      { 
        Px = px;
        Py = py;
        Pz = pz; 
      };
      inline G4ThreeVector GetMomentum() const
      { return G4ThreeVector(Px,Py,Pz); };
      inline G4double GetPx() const
      { return Px; };
      inline G4double GetPy() const
      { return Py; };
      inline G4double GetPz() const
      { return Pz; };
      inline void SetNext(G4PrimaryParticle * np)
      { 
        if(nextParticle == NULL)
        { nextParticle = np; }
        else
        { nextParticle->SetNext(np); }
      };
      inline G4PrimaryParticle * GetNext() const
      { return nextParticle; };
      inline void SetDaughter(G4PrimaryParticle * np)
      { 
        if(daughterParticle == NULL)
        { daughterParticle = np; }
        else
        { daughterParticle->SetNext(np); }
      };
      inline G4PrimaryParticle * GetDaughter() const
      { return daughterParticle; };
      inline void SetTrackID(G4int id)
      { trackID = id; };
      inline G4int GetTrackID() const
      { return trackID; };
      inline void SetMass(G4double mas)
      { mass = mas; };
      inline G4double GetMass() const
      { return mass; };
      inline void SetPolarization(G4double px,G4double py,G4double pz)
      {
        polX = px;
        polY = py;
        polZ = pz;
      };
      inline G4ThreeVector GetPolarization() const
      { return G4ThreeVector(polX,polY,polZ); };
      inline G4double GetPolX() const { return polX; };
      inline G4double GetPolY() const { return polY; };
      inline G4double GetPolZ() const { return polZ; };
};

extern G4Allocator<G4PrimaryParticle> aPrimaryParticleAllocator;

inline void * G4PrimaryParticle::operator new(size_t)
{
  void * aPrimaryParticle;
  aPrimaryParticle = (void *) aPrimaryParticleAllocator.MallocSingle();
  return aPrimaryParticle;
}

inline void G4PrimaryParticle::operator delete(void * aPrimaryParticle)
{
  aPrimaryParticleAllocator.FreeSingle((G4PrimaryParticle *) aPrimaryParticle);
}


#endif

