#ifndef GhadTrack_h
#define GhadTrack_h

#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4KineticTrack.hh"

class GhadTrack 
{
  public:
    GhadTrack(G4KineticTrack * anO) : org(anO) 
    {
      theType = anO->GetDefinition();
      theMom = anO->Get4Momentum();
      thePosition = anO->GetPosition();
    }
    
    G4bool operator == (const GhadTrack * aT)
    {
      if(this == aT) return true;
      return false;
    }
    
    // linear transport for time aStep
    void Go(G4double aStep)
    {
      // v = pc2/e
      G4cout << "Starting with "<<theMom.t()<<endl;
      G4cout << "in direction "<<theMom.vect().unit()<<endl;
      G4ThreeVector vel = c_light*theMom.vect()/theMom.t();
      G4double minStep = 0.5*fermi;
      G4ThreeVector d = vel*aStep;
      if(d.mag()<minStep)  d=vel*minStep/c_light;
      G4cout << "Old position "<<thePosition<<" "<<d<<endl;
      thePosition += d;
      G4cout << "Goint by "<<vel<<" "<<d<<" "<<d.mag()<<endl;
      G4cout << "New position "<<thePosition<<endl;
    }
    
    G4ThreeVector GetVelocity()
    {
      return c_light*theMom.vect()/theMom.t();
    }
    
    G4ThreeVector GetPosition() {return thePosition;}
    void SetPosition(G4ThreeVector & aV)
    {
      thePosition = aV;
    }
    
    G4KineticTrack * GetOrg() {return org;}
    G4ParticleDefinition * GetDefinition() {return theType;}
    G4LorentzVector & GetMom() {return theMom;}
    G4double GetActualMass() {return org->GetActualMass();}
    
    
  private:
  
  G4ParticleDefinition * theType;
  G4LorentzVector theMom;
  G4ThreeVector thePosition;
  G4KineticTrack * org;
};

#endif
