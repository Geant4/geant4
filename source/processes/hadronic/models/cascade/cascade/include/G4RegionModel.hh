//This class represents one region in a nucleus of the Bertini INC model.
//Hides the properties of the region and the generation of interaction partner. 
//The radiuses and densities are not the original.

#ifndef G4REGIONMODEL
#define G4REGIONMODEL

#include "G4ios.hh"
#include <vector>
#include <math.h>
#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4VRegionModel.hh"
#include "Randomize.hh"
typedef std::vector<G4double>::const_iterator iterator; 

class G4RegionModel : public G4VRegionModel
{
public:
  
  G4RegionModel();
  ~G4RegionModel();

  void CreateRegion(G4int A, G4int Z);
  void PrintRegion();

  G4double GetOuterRadius(); 
  G4double GetProtonDensity();
  G4double GetNeutronDensity();
  G4double GetProtonPotentialEnergy();
  G4double GetNeutronPotentialEnergy();
  G4double GetProtonMaximumMomentum();
  G4double GetNeutronMaximumMomentum();
  G4double GetProtonMaximumEnergy();
  G4double GetNeutronMaximumEnergy();
  
  G4DynamicParticle* GenerateProton();
  G4DynamicParticle* GenerateNeutron();    

private:
  G4int nucleusA; //A of the whole nucleus
  G4int nucleusZ;

  static G4int numberOfRegions;
  G4int myRegionNumber;
 
  G4double innerRadius;
  G4double outerRadius;
  G4double nuclearRadius;

  G4double protonDensity;
  G4double neutronDensity;

  G4double protonFermiMomentum;
  G4double neutronFermiMomentum;
  G4double protonFermiEnergy;
  G4double neutronFermiEnergy;

  G4double protonPotentialEnergy;
  G4double neutronPotentialEnergy;


  //methods:
  void InitializeRadius();
  void InitializeDensity();
  void InitializeFermi();
  void InitializePotentialEnergy();
  G4double GetFermiMomentum(G4double density, G4double mass);

};  
#endif








