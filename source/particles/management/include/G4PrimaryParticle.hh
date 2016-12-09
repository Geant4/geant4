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
// $Id: G4PrimaryParticle.hh 99159 2016-09-07 08:11:50Z gcosmo $
//
//


#ifndef G4PrimaryParticle_h
#define G4PrimaryParticle_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

#include "pwdefs.hh"

class G4ParticleDefinition;
class G4VUserPrimaryParticleInformation;

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
  
 public: // with description
  G4PrimaryParticle();
  G4PrimaryParticle(G4int Pcode);
  G4PrimaryParticle(G4int Pcode,
		    G4double px,G4double py,G4double pz);
  G4PrimaryParticle(G4int Pcode,
		    G4double px,G4double py,G4double pz,G4double E);
  G4PrimaryParticle(const G4ParticleDefinition* Gcode);
  G4PrimaryParticle(const G4ParticleDefinition* Gcode,
		    G4double px,G4double py,G4double pz);
  G4PrimaryParticle(const G4ParticleDefinition* Gcode,
		    G4double px,G4double py,G4double pz,G4double E);
  
 public:
  virtual ~G4PrimaryParticle();
  
  // copy constructor and assignment operator 
  // nextParticle and daughterParticle is copied by object (i.e. deep copy)
  // userInfo will nt be copied
  G4PrimaryParticle(const G4PrimaryParticle &right);   
  G4PrimaryParticle & operator=(const G4PrimaryParticle &right);
  
  // equal operator  returns 'true' only if the same object (i.e comarison by pointer value)
  G4int operator==(const G4PrimaryParticle &right) const;
  G4int operator!=(const G4PrimaryParticle &right) const;
  
 public: // with description
  void Print() const;
  // Print the properties of the particle.
  

 public: // with description
  // followings are get/set methods available.
  //   "trackID" will be set if this particle is sent to G4EventManager 
  //    and converted to G4Track. Otherwise = -1.
  //    The mass and charge in G4ParticleDefinition will be used in default.
  //   "SetMass" and "SetCharge" methods are used to set dynamical mass and charge 
  //   G4DynamicParticle."GetMass" and "GetCharge" methods will return 
  //   those in G4ParticleDefinition if users do not set dynamical ones. 

  G4int GetPDGcode()  const;
  void SetPDGcode(G4int Pcode);
  G4ParticleDefinition * GetG4code() const;
  void SetG4code(const G4ParticleDefinition * Gcode);
  const G4ParticleDefinition * GetParticleDefinition() const;
  void SetParticleDefinition(const G4ParticleDefinition * pdef);
  G4double GetMass() const;
  void SetMass(G4double mas);
  G4double GetCharge() const;
  void SetCharge(G4double chg);
  G4double GetKineticEnergy() const;
  void SetKineticEnergy(G4double eKin); 
  const G4ThreeVector& GetMomentumDirection() const;
  void SetMomentumDirection(const G4ThreeVector& p); 
  G4double GetTotalMomentum() const;
  void Set4Momentum(G4double px, G4double py, G4double pz, G4double E);
  G4double GetTotalEnergy() const;
  void SetTotalEnergy(G4double eTot ); 
  G4ThreeVector GetMomentum() const;
  void SetMomentum(G4double px, G4double py, G4double pz);
  G4double GetPx() const;
  G4double GetPy() const;
  G4double GetPz() const;
  G4PrimaryParticle * GetNext() const;
  void SetNext(G4PrimaryParticle * np);
  void ClearNext();
  G4PrimaryParticle * GetDaughter() const;
  void SetDaughter(G4PrimaryParticle * np);
  G4int GetTrackID() const;
  void SetTrackID(G4int id);
  G4ThreeVector GetPolarization() const;
  void SetPolarization(const G4ThreeVector& pol);
  void SetPolarization(G4double px,G4double py,G4double pz);
  G4double GetPolX() const;
  G4double GetPolY() const;
  G4double GetPolZ() const;
  G4double GetWeight() const;
  void SetWeight(G4double w);
  G4double GetProperTime() const;
  void SetProperTime(G4double t);
  G4VUserPrimaryParticleInformation* GetUserInformation() const;
  void SetUserInformation(G4VUserPrimaryParticleInformation* anInfo);

 private:
  G4int PDGcode;
  const G4ParticleDefinition * G4code;
  
  G4ThreeVector direction;
  G4double kinE;
  
  G4PrimaryParticle * nextParticle;
  G4PrimaryParticle * daughterParticle;
  
  G4int trackID;  // This will be set if this particle is
  // sent to G4EventManager and converted to
  // G4Track. Otherwise = -1.
  
  G4double mass;  
  G4double charge;
  G4double polX;
  G4double polY;
  G4double polZ;
  G4double Weight0;
  G4double properTime;
  G4VUserPrimaryParticleInformation* userInfo;

};

extern G4PART_DLL G4ThreadLocal G4Allocator<G4PrimaryParticle> *aPrimaryParticleAllocator;

inline void * G4PrimaryParticle::operator new(size_t)
{
  if (!aPrimaryParticleAllocator)
  {
    aPrimaryParticleAllocator = new G4Allocator<G4PrimaryParticle>;
  }
  return (void *) aPrimaryParticleAllocator->MallocSingle();
}

inline void G4PrimaryParticle::operator delete(void * aPrimaryParticle)
{
  aPrimaryParticleAllocator->FreeSingle((G4PrimaryParticle *) aPrimaryParticle);
}

inline G4double G4PrimaryParticle::GetMass() const
{ return mass;  }

inline G4double G4PrimaryParticle::GetCharge() const
{ return charge; }

inline G4int G4PrimaryParticle::GetPDGcode() const
{ return PDGcode; }
     
inline G4ParticleDefinition * G4PrimaryParticle::GetG4code() const
{ return const_cast<G4ParticleDefinition*>(G4code); }

inline const G4ParticleDefinition * G4PrimaryParticle::GetParticleDefinition() const
{ return G4code; }
    
inline G4double G4PrimaryParticle::GetTotalMomentum() const
{ 
  if (mass<0.)  return kinE; 
  else          return std::sqrt(kinE*(kinE+2.*mass));
}

inline G4ThreeVector G4PrimaryParticle::GetMomentum() const
{ return GetTotalMomentum()*direction;}

inline const G4ThreeVector& G4PrimaryParticle::GetMomentumDirection() const
{ return direction;}

inline void G4PrimaryParticle::SetMomentumDirection(const G4ThreeVector& p) 
{ direction = p;}

inline G4double G4PrimaryParticle::GetPx() const
{ return GetTotalMomentum()*direction.x();}

inline G4double G4PrimaryParticle::GetPy() const
{ return GetTotalMomentum()*direction.y();}

inline G4double G4PrimaryParticle::GetPz() const
{ return GetTotalMomentum()*direction.z();}

inline G4double G4PrimaryParticle::GetTotalEnergy() const
{ 
  if (mass<0.)  return kinE; 
  else          return kinE+mass;
}

inline void G4PrimaryParticle::SetTotalEnergy(G4double eTot ) 
{ 
  if (mass<0.)  kinE = eTot; 
  else          kinE = eTot - mass;
}

inline G4double G4PrimaryParticle::GetKineticEnergy() const
{ return kinE; }
   
inline void G4PrimaryParticle::SetKineticEnergy(G4double eKin) 
{ kinE = eKin; }

inline G4PrimaryParticle * G4PrimaryParticle::GetNext() const
{ return nextParticle; }

inline G4PrimaryParticle * G4PrimaryParticle::GetDaughter() const
{ return daughterParticle; }

inline G4int G4PrimaryParticle::GetTrackID() const
{ return trackID; }

inline G4ThreeVector G4PrimaryParticle::GetPolarization() const
{ return G4ThreeVector(polX,polY,polZ); }

inline G4double G4PrimaryParticle::GetPolX() const 
{ return polX; }

inline G4double G4PrimaryParticle::GetPolY() const 
{ return polY; }

inline G4double G4PrimaryParticle::GetPolZ() const 
{ return polZ; }

inline G4double G4PrimaryParticle::GetWeight() const 
{ return Weight0; }

inline void G4PrimaryParticle::SetWeight(G4double w) 
{ Weight0 = w; }

inline void G4PrimaryParticle::SetProperTime(G4double t)
{ properTime = t; }

inline G4double G4PrimaryParticle::GetProperTime() const 
{ return properTime; }

inline void G4PrimaryParticle::SetUserInformation(G4VUserPrimaryParticleInformation* anInfo)
{ userInfo = anInfo; }

inline G4VUserPrimaryParticleInformation* G4PrimaryParticle::GetUserInformation() const
{ return userInfo; }

inline void G4PrimaryParticle::SetG4code(const G4ParticleDefinition* Gcode)
{
  SetParticleDefinition(Gcode);
}

inline void G4PrimaryParticle::SetNext(G4PrimaryParticle * np)
{ 
  if   (nextParticle == 0) { nextParticle = np; }
  else                     { nextParticle->SetNext(np); }
}

inline void G4PrimaryParticle::ClearNext()
{
  nextParticle = nullptr;
}

inline void G4PrimaryParticle::SetDaughter(G4PrimaryParticle * np)
{ 
  if(daughterParticle == 0)  { daughterParticle = np; }
  else                       { daughterParticle->SetNext(np); }
}
     
inline void G4PrimaryParticle::SetTrackID(G4int id)
{ trackID = id; }

inline void G4PrimaryParticle::SetMass(G4double mas)
{ mass = mas; }

inline void G4PrimaryParticle::SetCharge(G4double chg)
{ charge = chg; }

inline void G4PrimaryParticle::SetPolarization(G4double px,G4double py,G4double pz)
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
