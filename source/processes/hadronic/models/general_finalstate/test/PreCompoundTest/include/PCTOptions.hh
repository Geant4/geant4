#ifndef PCTOptions_hh
#define PCTOptions_hh 1

#include "globals.hh"

#include <vector>

#include "PCTProjectileDirection.hh"
#include "PCTModes.hh"

class PCTOptions
{
public:
  inline PCTOptions();
  inline ~PCTOptions();
  
  void Initialize();
  
  inline G4String GetFilename() const;
  inline PCTEVAPModes GetEvapMode() const;
  inline PCTPREEQEmissionModes GetPreeqEmissionMode() const;
  inline PCTPREEQTransitionModes GetPreeqTransitionMode() const;
  inline PCTProjectileDirection GetProjectileDirection() const;
  inline G4int GetProjectileA() const;
  inline G4int GetProjectileZ() const;
  inline G4double GetProjectileKineticEnergy() const;
  inline G4int GetTargetA() const;
  inline G4int GetTargetZ() const;
  inline G4int GetNumberOfIsotopes() const;
  inline std::pair<G4int,G4double> GetIsotope(const unsigned int i) const;
  inline const std::vector< std::pair<G4int,G4double> >& GetTargetMaterial() const;
  inline G4bool IsNaturalTarget() const;
  inline G4int GetNumberOfIterations() const;
  inline G4int GetNumberOfParticles() const;
  inline G4int GetNumberOfHoles() const;
  inline G4int GetNumberOfCharged() const;
  inline const G4bool UsingINC() const;
  inline const G4bool UsingFermi() const;
  
private:
  inline PCTOptions(const PCTOptions& right);
  inline PCTOptions& operator=(const PCTOptions& right);
  
  inline void PutEmptyLine(const G4int n = 1);

  void AskForFermiBreakUpMode();
  void AskForPreeqEmissionMode();
  void AskForPreeqTransitionMode();
  void AskForEvaporationMode();
  void AskForTarget();
  void AskForProjectile();
  void AskForINC();
  void AskForExcitons();
  void AskForIterations();
  void GenerateFilename();
  
private:
  PCTEVAPModes theEvapMode;
  PCTPREEQEmissionModes thePreeqEmissionMode;
  PCTPREEQTransitionModes thePreeqTransitionMode;
  PCTProjectileDirection pdirection;
  G4int pA;
  G4int pZ;
  G4double pKineticEnergy;
  G4bool natural;
  G4int tA;
  G4int tZ;
  std::vector< std::pair<G4int,G4double> > material; 
  G4int iterations;
  G4int particles;
  G4int holes;
  G4int charged;
  G4bool inc;
  G4bool fermi;
  G4String filename;
};

#include "PCTOptions.icc"

#endif  // PCTOptions_hh
