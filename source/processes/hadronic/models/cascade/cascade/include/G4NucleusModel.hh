#ifndef G4NUCLEUSMODEL
#define G4NUCLEUSMODEL

#include <iostream.h>
#include "G4ThreeVector.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4RegionModel.hh"
#include "globals.hh"
#include "Randomize.hh"

enum G4ParticleType { G4Proton=0, G4Neutron=1 };
enum { NUMBER_OF_REGIONS = 1 };

class G4NucleusModel {
public:


  G4NucleusModel(const G4int A, const G4int Z);
  ~G4NucleusModel();

  void SetParameters(const G4int A, const G4int Z);

  G4double GetMassNumber();    
  G4double GetProtonNumber();
  G4double GetNuclearRadius(); 
  G4double GetNuclearMass();
  G4double GetAtomicMass();
  G4double GetExcitationEnergy();
  G4double GetNuclearDensity(G4ThreeVector position); //uses region model
  void AddExcitationEnergy(G4double additionalEnergy);

   //initializes incident particle's properties
  void DoCollision(G4ThreeVector initialMomentum, G4ParticleType particle); 
  G4ThreeVector GetInterActionPoint();
  G4bool IsInside(G4ThreeVector position);
  G4bool PauliPrinciple(G4double energy, G4double mass); //check Pauli blocking
  G4ParticleType ReturnTargetNucleon();//
  G4ThreeVector GetTargetMomentum();//
 
private:

//variables:
  G4double massNumber;
  G4double protonNumber;
  G4double radius;
  static const G4double radius0 = 1.0E-15;

  G4double bindingEnergy;
  G4double  excitationEnergy;   
  G4RegionModel* shellStructure;

  //impact variables:
  G4ThreeVector incidencePoint;
  G4double interActionLength;
  G4ThreeVector interActionPoint;
  G4ParticleType targetParticle;
  G4ThreeVector targetMomentum;
  G4ParticleType incidentParticle;
  G4ThreeVector incidentMomentum;
  G4double incidentEnergy;

//methods:
  G4ThreeVector GetIncidencePoint();
  G4double GetInterActionLength();
  G4bool IsBoundaryCrossed();
  G4ThreeVector GetCrossingPoint();
  //changes the direction of incident particle, if it crosses a region boundary
  void doReflection(); 
};
#endif



