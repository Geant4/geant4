#ifndef G4HadSecondary_hh
#define G4HadSecondary_hh

#include "G4DynamicParticle.hh"

class G4HadSecondary
{
  public:
    G4HadSecondary(G4DynamicParticle * aT, G4double aWeight = 1);
    G4DynamicParticle * GetParticle() {return theP;}
    G4double GetWeight() {return theWeight;}
    void SetWeight(G4double aW){theWeight= aW;}
    void SetTime(G4double aT) {theTime = aT;}
    G4double GetTime() {return theTime;}
    
    
  private:
  
   G4HadSecondary(){};
   
   G4DynamicParticle * theP; 
   G4double theWeight;
   G4double theTime;
};

#endif
