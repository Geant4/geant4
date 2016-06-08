// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PPrimaryParticle.ddl,v 1.2 2000/11/02 12:41:59 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

#ifndef G4PPrimaryParticle_h
#define G4PPrimaryParticle_h 1

#include "G4Pglobals.hh"
#include "G4ThreeVector.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PrimaryParticle;

// class description:
//
// This is the persistent version of a class G4PrimaryParticle.
//

class G4PPrimaryParticle 
 : public HepPersObj
{
  public:
      G4PPrimaryParticle( G4PrimaryParticle* particle );
      ~G4PPrimaryParticle();

      G4PrimaryParticle*  MakeTransientObject();

  private:
      G4Pint PDGcode;
      // G4ParticleDefinition * G4code;
      G4Pdouble Px;
      G4Pdouble Py;
      G4Pdouble Pz;
      
      d_Ref<G4PPrimaryParticle> nextParticle;
      d_Ref<G4PPrimaryParticle> daughterParticle;

      G4Pint trackID;  // This will be set if this particle is
                       // sent to G4EventManager and converted to
                       // G4Track. Otherwise = -1.
      G4Pdouble mass;  // This is just for book keeping.
                       // This will not be used but the mass in
                       // G4ParticleDefinition will be used.

      G4Pdouble polX;
      G4Pdouble polY;
      G4Pdouble polZ;

  public: // with description
      // followings are get methods available.
      //   "trackID" will be set if this particle is sent to G4EventManager and converted to
      //   G4Track. Otherwise = -1.
      //   "mass" is just for book keeping. This will not be used but the mass in
      //   G4ParticleDefinition will be used.
      inline G4int GetPDGcode() const
      { return PDGcode; }
      // inline G4ParticleDefinition * GetG4code() const
      // { return G4code; }
      inline G4ThreeVector GetMomentum() const
      { return G4ThreeVector(Px,Py,Pz); }
      inline G4double GetPx() const
      { return Px; }
      inline G4double GetPy() const
      { return Py; }
      inline G4double GetPz() const
      { return Pz; }
      inline HepRef(G4PPrimaryParticle) GetNext() const
      { return nextParticle; }
      inline HepRef(G4PPrimaryParticle) GetDaughter() const
      { return daughterParticle; }
      inline G4int GetTrackID() const
      { return trackID; }
      inline G4double GetMass() const
      { return mass; }
      inline G4ThreeVector GetPolarization() const
      { return G4ThreeVector(polX,polY,polZ); }
      inline G4double GetPolX() const { return polX; }
      inline G4double GetPolY() const { return polY; }
      inline G4double GetPolZ() const { return polZ; }

};

#endif

