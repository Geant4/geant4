#ifndef G4REGIONMODEL
#define G4REGIONMODEL

#include "G4ios.hh"
#include <vector>
#include <math.h>
#include "globals.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
//#include "G4NucleusModel.hh"


typedef std::vector<G4double>::const_iterator iterator; 

class G4RegionModel //:public G4VRegionModel
{
public:
  G4RegionModel(const G4int numberOfLayers, const G4int A, const G4int Z);
  ~G4RegionModel();

  //instead of A and Z outer radius of the nucleus?
  //void Init(const G4int numberOfLayers, const G4int A, const G4int Z);
  G4double GetDensity(G4double radius);
  G4double GetPotentialEnergy(G4double r, G4int particle);
  G4double GetMaximumNucleonMomentum(G4double radius, G4int nucleon);
  // G4double NumberOfRegions();   

private:

  G4int massNumber;
  G4int protonNumber;
  vector<G4double> radius; //contains the outer radiuses of the shells
  vector<G4double> density;

  vector<G4double> protonFermiEnergy;
  vector<G4double> neutronFermiEnergy;
  vector<G4double> protonFermiMomentum;
  vector<G4double> neutronFermiMomentum;
  
  vector<G4double> protonPotentialEnergy;
  vector<G4double> neutronPotentialEnergy;

  static const G4double radius0=1.0E-15; 
  static const G4double BE = 7;
  //static const G4double pi = 3.141592;

  G4double GetFermiMomentum(G4double density, G4double mass);
  G4double GetFermiEnergy(G4double density, G4double mass);
};  
#endif








