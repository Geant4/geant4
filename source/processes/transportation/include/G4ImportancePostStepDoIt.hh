#ifndef G4ImportancePostStepDoIt_hh
#define G4ImportancePostStepDoIt_hh G4ImportancePostStepDoIt_hh

class G4VImportanceSampler;
class G4ParticleChange;
class G4Track;
class G4Step;
class G4Nsplit_Weight;

class G4ImportancePostStepDoIt{
public:
  G4ImportancePostStepDoIt(){}
  ~G4ImportancePostStepDoIt(){}
  
  void DoIt(const G4Track& aTrack, 
	    G4ParticleChange *aParticleChange, 
	    const G4Nsplit_Weight nw);
  
private:
  void Split(const G4Track &aTrack,
	     const G4Nsplit_Weight &nw,
	     G4ParticleChange *aParticleChange);
};



#endif
