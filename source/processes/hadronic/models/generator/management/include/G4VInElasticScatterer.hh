#ifndef G4VInElasticScatterer_h
#define G4VInElasticScatterer_h 1

class G4KineticTrackVector;
class G4KineticTrack;

class G4VInElasticScatterer 
   {
public:

   G4VInElasticScatterer();
   virtual ~G4VInElasticScatterer();
    
public:    
    virtual G4KineticTrackVector* Scatter(const G4KineticTrack &aProjectile, const G4KineticTrack &aTarget) = 0;

};
  
#endif 


