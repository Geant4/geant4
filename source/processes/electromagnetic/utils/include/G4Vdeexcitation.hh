#ifndef G4VDEEXCITATION
#define G4VDEEXCITATION 1

#include "globals.hh"

using namespace std;

class G4VDeexcitation {

  //initialization
  virtual void PreparePhysicsTable();
  virtual void BuildPhysicsTable();

  // this is needed for the Lowenergy processes (photoelectric and eIoni)
  virtual G4AtomicChell* GetAtomicShell(G4int Z, G4int ShellIndex);


  virtual void AlongStepDeexcitation(vector<G4DynamicParticle*>* secVect, 
				     G4DynamicParticle* icidentParticle, G4MaterialCutsCouple, 
				     G4double trueStepLenght, G4double eLoss);

  // method for having deexcitation producta
  virtual void GenerateParticles(vector<G4DynamicParticle*>*, G4int Z, G4double);

  virtual G4bool CheckActiveRegion(G4int coupleIndex);
  virtual G4bool IsFluorescenceActive(); // maybe we need to specify region
  virtual G4bool IsPIXECrossSectionActive(); // maybe we need to specify region

  virtual void SetFluorescenceActiveRegion(G4Region region=0);
  virtual void SetAugerActiveRegion(G4Region region=0);
  // need to know if CS must be caluclated for this region
  virtual void SetPIXECrossSectionActiveRegion(G4Region region=0); 

  virtual void SetPIXECrossSectionModel(G4String name);

  // method used to select random sehll according to models.
  virtual G4AtomicShell* SelectRandomShell(G4int Z, G4DynamicParticle*);

  // method used to get PIXE CS in the tables created at initialization time
  virtual G4double GetPIXECrossSection (G4int Z, G4ParticleDefinition*, G4double kinE);
  // method used to calculate PIXE CS directly from models
  virtual G4double CalculatePIXECrossSection (G4int Z, G4ParticleDefinition*, G4double kinE);

};

#endif

