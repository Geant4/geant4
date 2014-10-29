#ifndef G4LEPTSElasticModel_h
#define G4LEPTSElasticModel_h 1

#include "G4VLEPTSModel.hh"
#include "G4ParticleChangeForGamma.hh"

class G4LEPTSElasticModel : public G4VLEPTSModel { 
public:
  G4LEPTSElasticModel(const G4String& modelName ="G4LEPTSElasticModel");
  ~G4LEPTSElasticModel();

  virtual void Initialise(const G4ParticleDefinition*, 
                          const G4DataVector&);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin = 0.0,
                                 G4double tmax = DBL_MAX);

 // main method to compute cross section per Volume
  virtual G4double CrossSectionPerVolume(const G4Material*,
                                         const G4ParticleDefinition*,
                                         G4double kineticEnergy,
                                         G4double cutEnergy = 0.0,
                                         G4double maxEnergy = DBL_MAX);

private:

  std::map<const G4Material*, G4double> theMassTarget;//M*c2 target, proj
  std::map<const G4Material*, G4double> theMassProjectile; 
  G4double EnergyTransfer(G4double, G4double, G4double, G4double);

private:
  G4ParticleChangeForGamma* fParticleChangeForGamma;
};


#endif
