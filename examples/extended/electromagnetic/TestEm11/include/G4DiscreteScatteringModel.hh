#ifndef G4DiscreteScatteringModel_HH
#define G4DiscreteScatteringModel_HH

#include "G4VEmModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementData.hh"
#include "G4PhysicsVector.hh"
#include "G4ParticleChangeForGamma.hh"

#include <vector>

typedef std::vector<G4double> G4PVDataVector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4DiscreteScatteringModel : public G4VEmModel
{

public:
  
  G4DiscreteScatteringModel(G4int iNumAngles=1);

  virtual ~G4DiscreteScatteringModel();
  
  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  
  void ReadData(G4int, const G4String & argFileName);
          
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                G4double E, 
                G4double Z, 
                G4double A = 0., 
                G4double cutEnergy = 0.0,
                G4double maxEnergy = DBL_MAX);
                  
  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
                                     const G4MaterialCutsCouple*, 
                                     const G4DynamicParticle*, 
                                     G4double tmin=0.0, 
                                     G4double maxEnergy=DBL_MAX);
  
  inline void SetNumberOfAngles(G4int N)        { fNumAngles=N; };
  inline void SetAnalog(const G4String& model)  { fAnalogModel=model; };
      
private:

  G4ThreeVector GetNewDirection(G4double z1);

  static G4ElementData*     fCdf;
  static G4ElementData*     fTcs;
  G4PVDataVector            fGrid;
  G4ParticleChangeForGamma* fParticleChange; 
  G4String                  fAnalogModel;
  G4int                     fNumAngles;
  G4double                  fLowEnergyLimit;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


