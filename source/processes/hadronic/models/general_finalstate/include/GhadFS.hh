#ifndef GhadFS_h
#define GhadFS_h

enum GhadAction {GEO=1, HIT=2, FLIGHT=0};

#include <vector>
#include "GhadTrack.hh"
#include "G4Nucleon.hh"
#include "GhadNucleus.hh"


class GhadFS : public std::vector<GhadTrack>
{
  public:
    GhadFS(GhadAction & act, G4double aStep, GhadParticles::iterator& aPro, G4Nucleon * aT) 
          : didAct(act), theStep(aStep), thePro(aPro), theTarget(aT) {}
	  
    virtual ~GhadFS() {}
    
    G4double GetTime() { return theStep; }
    
    G4bool ImVoid() { return didAct==FLIGHT; }
    
    GhadParticles::iterator & GetProjectile()
    {
      if(didAct) return thePro;
      else G4Exception("GhadFS: getting projectile while nothing had happened");
    }
    
    void HitTarget() {theTarget->Hit();}
  
  private:
    GhadAction didAct;
    G4double theStep;
    GhadParticles::iterator & thePro;
    G4Nucleon * theTarget;
};

#endif
