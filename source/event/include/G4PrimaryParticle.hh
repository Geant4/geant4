// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PrimaryParticle.hh,v 1.5 2000-10-18 12:41:21 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//


#ifndef G4PrimaryParticle_h
#define G4PrimaryParticle_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class G4ParticleDefinition;

// class description:
//
//  This is a class which represents a primary particle.
// This is completely deferent class from G4Track or G4DynamicParticle.
// This class is designed with taking into account the possibility of
// making this class persistent, i.e. kept with G4Event class object
// to ODBMS. Thus this class is almost free from any other Geant4 classes.
// The only exception is a pointer to G4ParticleDefinition but it can be
// rebuilt by the PDGcode.
//
//  Primary particles are stored in G4PrimaryVertex object with a form
// of linked list. Also, an object of this PrimaryParticle class can have
// one or more objects of this class as its daughters with a form of 
// linked list.
//  A parimary particle represented by this class object needs not to be
// a particle of type which Geant4 can simulate.
//  case a) mother particle is not a particle Geant4 can simulate
//   daughters associated to the mother will be examined.
//  case b) mother particle is a perticle Geant4 can simulate
//   daughters associated to the mother will be converted to G4Dynamic 
//   particle and be set to the mother G4Track. For this case, dauthers
//   are used as the "pre-fixed" decay channel.
//

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

  public: // with description
      void Print() const;
      // Print the properties of the particle.

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
      G4double charge;
      G4double polX;
      G4double polY;
      G4double polZ;

  public: // with description
      // followings are get methods available.
      //   "trackID" will be set if this particle is sent to G4EventManager and converted to
      //   G4Track. Otherwise = -1.
      //   "mass" is just for book keeping. This will not be used but the mass in
      //   G4ParticleDefinition will be used.
      inline G4int GetPDGcode() const
      { return PDGcode; }
      inline G4ParticleDefinition * GetG4code() const
      { return G4code; }
      inline G4ThreeVector GetMomentum() const
      { return G4ThreeVector(Px,Py,Pz); }
      inline G4double GetPx() const
      { return Px; }
      inline G4double GetPy() const
      { return Py; }
      inline G4double GetPz() const
      { return Pz; }
      inline G4PrimaryParticle * GetNext() const
      { return nextParticle; }
      inline G4PrimaryParticle * GetDaughter() const
      { return daughterParticle; }
      inline G4int GetTrackID() const
      { return trackID; }
      inline G4double GetMass() const
      { return mass; }
      inline G4double GetCharge() const
      { return charge; }
      inline G4ThreeVector GetPolarization() const
      { return G4ThreeVector(polX,polY,polZ); }
      inline G4double GetPolX() const { return polX; }
      inline G4double GetPolY() const { return polY; }
      inline G4double GetPolZ() const { return polZ; }

  public: // with description
      // Followings are available Set methods.
      void SetPDGcode(G4int Pcode);
      void SetG4code(G4ParticleDefinition * Gcode);
      inline void SetMomentum(G4double px, G4double py, G4double pz)
      { 
        Px = px;
        Py = py;
        Pz = pz; 
      }
      inline void SetNext(G4PrimaryParticle * np)
      { 
        if(nextParticle == NULL)
        { nextParticle = np; }
        else
        { nextParticle->SetNext(np); }
      }
      inline void SetDaughter(G4PrimaryParticle * np)
      { 
        if(daughterParticle == NULL)
        { daughterParticle = np; }
        else
        { daughterParticle->SetNext(np); }
      }
      inline void SetTrackID(G4int id)
      { trackID = id; }
      inline void SetMass(G4double mas)
      { mass = mas; }
      inline void SetCharge(G4double chg)
      { charge = chg; }
     inline void SetPolarization(G4double px,G4double py,G4double pz)
      {
        polX = px;
        polY = py;
        polZ = pz;
      }
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

