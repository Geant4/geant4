#ifndef G4DNADissociation_h
#define G4DNADissociation_h 1

#include "G4VEmProcess.hh"

class G4DNADissociation : public G4VEmProcess
{ 
public:
  G4DNADissociation(const G4String& processName ="G4DNADissociation");
  ~G4DNADissociation();

  G4bool IsApplicable(const G4ParticleDefinition&);
  virtual void InitialiseProcess(const G4ParticleDefinition*);
  virtual void PrintInfo();
  
private:
  G4bool       isInitialised;

};


#endif
