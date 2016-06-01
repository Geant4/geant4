#ifndef G4VInElasticScatterer_h
#define G4VInElasticScatterer_h 1

#include "G4KineticTrackVector.hh"

class G4VInElasticScatterer 
   {
public:

   G4VInElasticScatterer();
   virtual ~G4VInElasticScatterer();
    
public:    
    virtual G4KineticTrackVector* Scatter(const G4KineticTrack &aProjectile, const G4KineticTrack &aTarget) = 0;

};
  
#endif 


