#ifndef HistoTest47_H
#define HistoTest47_H

#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"

#include <string>
#include <map>
#include <vector>

class G4VParticleChange;

class HistoTest47 {

public:

  HistoTest47(std::string namePart, std::string nameMat, G4double momentum,
	      std::string nameGen);
  virtual ~HistoTest47();

  void setParticle(std::string namePart);
  void setTarget(std::string nameMat)    {target = nameMat; unInitialized =true;}
  void setMomentum(double momentum)      {energy = momentum; unInitialized =true;}
  void setGenerator(std::string nameGen) {generator = nameGen; unInitialized =true;}
  void setJobID( int id ) { jobID = id ; return ; }
  virtual void fill(G4VParticleChange*, G4LorentzVector)=0;
  virtual void write(G4double cross_sec, G4int nevt)=0;

protected:

  G4int particleType(G4ParticleDefinition*);
  bool                                   unInitialized;
  std::string                            particle, target, generator;
  double                                 energy;
  int                                    jobID;

private:

  std::map<G4ParticleDefinition*, G4int> mapParticle;
  
};

#endif
