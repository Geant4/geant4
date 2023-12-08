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
// G4PrimaryParticle
//
// Class description:
//
// This class represents a primary particle.
// G4PrimaryParticle is a completely different object from G4Track or
// G4DynamicParticle. G4PrimaryParticle is designed to take into account
// the possibility of making the object persistent, i.e. kept with G4Event
// class object within an ODBMS. This class is therefore almost independent
// from any other Geant4 object. The only exception is a pointer to
// G4ParticleDefinition, which can be rebuilt by its PDGcode.
//
// Primary particles are stored in G4PrimaryVertex through a form of linked
// list. A G4PrimaryParticle object can have one or more objects of this
// class as its daughters as a form of linked list.
// A primary particle represented by this class object needs not to be
// a particle type which Geant4 can simulate:
//  Case a) mother particle is not a particle Geant4 can simulate
//          daughters associated to the mother will be examined.
//  Case b) mother particle is a perticle Geant4 can simulate
//          daughters associated to the mother will be converted to
//          G4DynamicParticle and be set to the mother G4Track.
//          For this case, daugthers are used as the "pre-fixed"
//          decay channel.

// Authors: G.Cosmo, 2 December 1995 - Design, based on object model
//          M.Asai, 29 January 1996 - First implementation
// --------------------------------------------------------------------
#ifndef G4PrimaryParticle_hh
#define G4PrimaryParticle_hh 1

#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "pwdefs.hh"

class G4ParticleDefinition;
class G4VUserPrimaryParticleInformation;

class G4PrimaryParticle
{
  public:
    // Constructors
    G4PrimaryParticle();
    G4PrimaryParticle(G4int Pcode);
    G4PrimaryParticle(G4int Pcode, G4double px, G4double py, G4double pz);
    G4PrimaryParticle(G4int Pcode, G4double px, G4double py, G4double pz, G4double E);
    G4PrimaryParticle(const G4ParticleDefinition* Gcode);
    G4PrimaryParticle(const G4ParticleDefinition* Gcode, G4double px, G4double py, G4double pz);
    G4PrimaryParticle(const G4ParticleDefinition* Gcode, G4double px, G4double py, G4double pz,
                      G4double E);

    // Destructor
    virtual ~G4PrimaryParticle();

    // Copy constructor and assignment operator
    // NOTE: nextParticle and daughterParticle are copied by object
    // (i.e. deep copy); userInfo will not be copied
    G4PrimaryParticle(const G4PrimaryParticle& right);
    G4PrimaryParticle& operator=(const G4PrimaryParticle& right);

    // Equality operator returns 'true' only if same object
    // (i.e. comparison by pointer value)
    G4bool operator==(const G4PrimaryParticle& right) const;
    G4bool operator!=(const G4PrimaryParticle& right) const;

    // Overloaded new/delete operators
    inline void* operator new(std::size_t);
    inline void operator delete(void* aPrimaryParticle);

    // Print the properties of the particle
    void Print() const;

    // Followings are the available accessors/modifiers.
    //   "trackID" will be set if this particle is sent to G4EventManager
    //    and converted to G4Track. Otherwise = -1.
    //    The mass and charge in G4ParticleDefinition will be used by default.
    //   "SetMass" and "SetCharge" methods are used to set dynamical mass and
    //    charge of G4DynamicParticle.
    //   "GetMass" and "GetCharge" methods will return those in
    //    G4ParticleDefinition if users do not set dynamical ones

    inline G4int GetPDGcode() const;
    void SetPDGcode(G4int Pcode);
    inline G4ParticleDefinition* GetG4code() const;
    inline void SetG4code(const G4ParticleDefinition* Gcode);
    inline const G4ParticleDefinition* GetParticleDefinition() const;
    void SetParticleDefinition(const G4ParticleDefinition* pdef);
    inline G4double GetMass() const;
    inline void SetMass(G4double mas);
    inline G4double GetCharge() const;
    inline void SetCharge(G4double chg);
    inline G4double GetKineticEnergy() const;
    inline void SetKineticEnergy(G4double eKin);
    inline const G4ThreeVector& GetMomentumDirection() const;
    inline void SetMomentumDirection(const G4ThreeVector& p);
    inline G4double GetTotalMomentum() const;
    void Set4Momentum(G4double px, G4double py, G4double pz, G4double E);
    inline G4double GetTotalEnergy() const;
    inline void SetTotalEnergy(G4double eTot);
    inline G4ThreeVector GetMomentum() const;
    void SetMomentum(G4double px, G4double py, G4double pz);
    inline G4double GetPx() const;
    inline G4double GetPy() const;
    inline G4double GetPz() const;
    inline G4PrimaryParticle* GetNext() const;
    inline void SetNext(G4PrimaryParticle* np);
    inline void ClearNext();
    inline G4PrimaryParticle* GetDaughter() const;
    inline void SetDaughter(G4PrimaryParticle* np);
    inline G4int GetTrackID() const;
    inline void SetTrackID(G4int id);
    inline G4ThreeVector GetPolarization() const;
    inline void SetPolarization(const G4ThreeVector& pol);
    inline void SetPolarization(G4double px, G4double py, G4double pz);
    inline G4double GetPolX() const;
    inline G4double GetPolY() const;
    inline G4double GetPolZ() const;
    inline G4double GetWeight() const;
    inline void SetWeight(G4double w);
    inline G4double GetProperTime() const;
    inline void SetProperTime(G4double t);
    inline G4VUserPrimaryParticleInformation* GetUserInformation() const;
    inline void SetUserInformation(G4VUserPrimaryParticleInformation* anInfo);

  private:
    const G4ParticleDefinition* G4code = nullptr;

    G4ThreeVector direction;
    G4double kinE = 0.0;

    G4PrimaryParticle* nextParticle = nullptr;
    G4PrimaryParticle* daughterParticle = nullptr;

    G4double mass = -1.0;
    G4double charge = 0.0;
    G4double polX = 0.0;
    G4double polY = 0.0;
    G4double polZ = 0.0;
    G4double Weight0 = 1.0;
    G4double properTime = -1.0;
    G4VUserPrimaryParticleInformation* userInfo = nullptr;

    G4int PDGcode = 0;
    G4int trackID = -1;  // This will be set if this particle is
                         // sent to G4EventManager and converted to
                         // G4Track. Otherwise = -1
};

extern G4PART_DLL G4Allocator<G4PrimaryParticle>*& aPrimaryParticleAllocator();

// ------------------------
// Inline methods
// ------------------------

inline void* G4PrimaryParticle::operator new(std::size_t)
{
  if (aPrimaryParticleAllocator() == nullptr) {
    aPrimaryParticleAllocator() = new G4Allocator<G4PrimaryParticle>;
  }
  return (void*)aPrimaryParticleAllocator()->MallocSingle();
}

inline void G4PrimaryParticle::operator delete(void* aPrimaryParticle)
{
  aPrimaryParticleAllocator()->FreeSingle((G4PrimaryParticle*)aPrimaryParticle);
}

inline G4double G4PrimaryParticle::GetMass() const
{
  return mass;
}

inline G4double G4PrimaryParticle::GetCharge() const
{
  return charge;
}

inline G4int G4PrimaryParticle::GetPDGcode() const
{
  return PDGcode;
}

inline G4ParticleDefinition* G4PrimaryParticle::GetG4code() const
{
  return const_cast<G4ParticleDefinition*>(G4code);
}

inline const G4ParticleDefinition* G4PrimaryParticle::GetParticleDefinition() const
{
  return G4code;
}

inline G4double G4PrimaryParticle::GetTotalMomentum() const
{
  if (mass < 0.) return kinE;
  return std::sqrt(kinE * (kinE + 2. * mass));
}

inline G4ThreeVector G4PrimaryParticle::GetMomentum() const
{
  return GetTotalMomentum() * direction;
}

inline const G4ThreeVector& G4PrimaryParticle::GetMomentumDirection() const
{
  return direction;
}

inline void G4PrimaryParticle::SetMomentumDirection(const G4ThreeVector& p)
{
  direction = p;
}

inline G4double G4PrimaryParticle::GetPx() const
{
  return GetTotalMomentum() * direction.x();
}

inline G4double G4PrimaryParticle::GetPy() const
{
  return GetTotalMomentum() * direction.y();
}

inline G4double G4PrimaryParticle::GetPz() const
{
  return GetTotalMomentum() * direction.z();
}

inline G4double G4PrimaryParticle::GetTotalEnergy() const
{
  if (mass < 0.) return kinE;
  return kinE + mass;
}

inline void G4PrimaryParticle::SetTotalEnergy(G4double eTot)
{
  if (mass < 0.)
    kinE = eTot;
  else
    kinE = eTot - mass;
}

inline G4double G4PrimaryParticle::GetKineticEnergy() const
{
  return kinE;
}

inline void G4PrimaryParticle::SetKineticEnergy(G4double eKin)
{
  kinE = eKin;
}

inline G4PrimaryParticle* G4PrimaryParticle::GetNext() const
{
  return nextParticle;
}

inline G4PrimaryParticle* G4PrimaryParticle::GetDaughter() const
{
  return daughterParticle;
}

inline G4int G4PrimaryParticle::GetTrackID() const
{
  return trackID;
}

inline G4ThreeVector G4PrimaryParticle::GetPolarization() const
{
  return G4ThreeVector(polX, polY, polZ);
}

inline G4double G4PrimaryParticle::GetPolX() const
{
  return polX;
}

inline G4double G4PrimaryParticle::GetPolY() const
{
  return polY;
}

inline G4double G4PrimaryParticle::GetPolZ() const
{
  return polZ;
}

inline G4double G4PrimaryParticle::GetWeight() const
{
  return Weight0;
}

inline void G4PrimaryParticle::SetWeight(G4double w)
{
  Weight0 = w;
}

inline void G4PrimaryParticle::SetProperTime(G4double t)
{
  properTime = t;
}

inline G4double G4PrimaryParticle::GetProperTime() const
{
  return properTime;
}

inline void G4PrimaryParticle::SetUserInformation(G4VUserPrimaryParticleInformation* anInfo)
{
  userInfo = anInfo;
}

inline G4VUserPrimaryParticleInformation* G4PrimaryParticle::GetUserInformation() const
{
  return userInfo;
}

inline void G4PrimaryParticle::SetG4code(const G4ParticleDefinition* Gcode)
{
  SetParticleDefinition(Gcode);
}

inline void G4PrimaryParticle::SetNext(G4PrimaryParticle* np)
{
  if (nextParticle == nullptr) {
    nextParticle = np;
  }
  else {
    nextParticle->SetNext(np);
  }
}

inline void G4PrimaryParticle::ClearNext()
{
  nextParticle = nullptr;
}

inline void G4PrimaryParticle::SetDaughter(G4PrimaryParticle* np)
{
  if (daughterParticle == nullptr) {
    daughterParticle = np;
  }
  else {
    daughterParticle->SetNext(np);
  }
}

inline void G4PrimaryParticle::SetTrackID(G4int id)
{
  trackID = id;
}

inline void G4PrimaryParticle::SetMass(G4double mas)
{
  mass = mas;
}

inline void G4PrimaryParticle::SetCharge(G4double chg)
{
  charge = chg;
}

inline void G4PrimaryParticle::SetPolarization(G4double px, G4double py, G4double pz)
{
  polX = px;
  polY = py;
  polZ = pz;
}

inline void G4PrimaryParticle::SetPolarization(const G4ThreeVector& pol)
{
  polX = pol.x();
  polY = pol.y();
  polZ = pol.z();
}

#endif
