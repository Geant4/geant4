//This class represents a nucleus and has all normal nuclear properties.
// It can also generate the interaction partner to the cascade particle.
//This class hides the structure of the nucleus and the calculation of impact variables from its clients. 
//The calculation of excitation energy is missing. 

#ifndef G4NUCLEUSMODEL
#define G4NUCLEUSMODEL

#include <iostream.h>
#include <string.h>
#include "G4ThreeVector.hh"
#include "G4NucleiPropertiesTable.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4Electron.hh"
#include "G4Deuteron.hh"
#include "G4Alpha.hh"
#include "G4DynamicParticle.hh"
#include "G4RegionModel.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4BertiniData.hh"
#include "G4LorentzVector.hh"
#include "G4KineticTrack.hh"

enum { NUMBER_OF_REGIONS = 3 }; //number of regions in the model

//definitions
typedef std::vector<G4VRegionModel*>::const_iterator regionIterator;
typedef std::vector<G4KineticTrack*>::const_iterator trackIterator;
typedef vector<G4KineticTrack*> trackVector;

class G4NucleusModel
{
public:

  G4NucleusModel();
  ~G4NucleusModel();
  void CreateModel(const G4int A, const G4int Z); //initializes the nucleus model

  void SetMomentum(G4ThreeVector momentum);
  void SetParameters(const G4int A, const G4int Z);
  void ChangeParameters(const G4int deltaA, const G4int deltaZ); 

  G4double GetA();    
  G4double GetZ();
  G4double GetRadius(); 
  G4double GetRegionRadius(G4int regionNumber);
  G4double GetNuclearMass();
  G4double GetAtomicMass();
  G4double GetBindingEnergy();
  G4double GetExcitationEnergy();
  G4double GetKineticEnergy();
  G4ThreeVector GetMomentum();

  G4double GetProtonDensity(G4ThreeVector position); 
  G4double GetNeutronDensity(G4ThreeVector position);
  G4double GetProtonFermiEnergy(G4ThreeVector position); //Fermi kinetic energy
  G4double GetNeutronFermiEnergy(G4ThreeVector position);  
  G4double GetProtonFermiMomentum(G4ThreeVector position);
  G4double GetNeutronFermiMomentum(G4ThreeVector position);
  G4double GetProtonPotentialEnergy(G4ThreeVector position);
  G4double GetNeutronPotentialEnergy(G4ThreeVector position); 
  
  void ChangeExcitationEnergy(G4double deltaE);
  void PrintModel(); 

   //calculates the recoil energy & momentum and takes them into account
  void ApplyRecoil(trackVector track);
 
  //initializes incident particle's properties and calculates impact parameters
  void DoCollision(G4DynamicParticle* particle); 
  void DoCollision(G4DynamicParticle* particle, G4ThreeVector coordinates);

  //can be used to do a boundary transition to an arbitrary particle,
  //changes particle's momentum
  G4VRegionModel* DoBoundaryTransition(G4VRegionModel* currentRegion,
				       G4ThreeVector crossingPoint,
				       G4DynamicParticle* particle);
 
  //calculates the point where 'particle' crosses two regions
  G4ThreeVector GetCrossingPoint(G4VRegionModel* currentRegion,
				 G4ThreeVector initialPoint,
				 G4DynamicParticle* particle);

  //returns a pointer to the region in radius 'r'
  G4VRegionModel* GetRegion(G4double r); 
  
  G4ThreeVector GetInteractionPoint();
  G4bool IsInside(G4ThreeVector position);

  //check Pauli blocking
  G4bool CheckPauliPrinciple(G4DynamicParticle *productParticle);

  //generates and returns the nucleon where hit
  G4DynamicParticle* ReturnTargetNucleon();
 
  //private:

//variables:
  G4double myA;
  G4double myZ;
  G4double radius;
  G4double atomicMass;

  G4double bindingEnergy;
  G4double excitationEnergy;
  G4double kineticEnergy;
  G4ThreeVector momentum; 
  vector<G4VRegionModel*> regions;
  static G4BertiniData* xSecDataBase;

//impact variables:
  G4ThreeVector incidencePoint;
  G4double interactionLength;
  G4ThreeVector interactionPoint;
  G4DynamicParticle* incidentParticle;
  G4bool movingIn;

//methods:
  G4double CalculateAtomicMass();
  G4ThreeVector GetIncidencePoint();
  G4double GetInteractionLength();
  G4ThreeVector CalculateInteractionPoint();

  //checks if a boundary is crossed
  G4bool BoundaryCrossed(G4VRegionModel* initialRegion,
			 G4ThreeVector initialPoint,	
			 G4DynamicParticle* particle,
			 G4double distance);

  //returns the pointer to next region where particle is heading for
  G4VRegionModel* GetNextRegion(G4VRegionModel* currentRegion,
			      G4ThreeVector initCoordinates,
			      G4DynamicParticle* particle);

  //excitation energy not working properly
  //G4double CalculateExcitationEnergy(); 
  //G4double normal();
  //void ApplyExcitationEnergy(); 
};
#endif






