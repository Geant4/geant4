#ifndef G4DNAPositronium_h
#define G4DNAPositronium_h 1

#include "G4VEmProcess.hh"

class G4DNAPositronium : public G4VEmProcess
{ 
public:
  G4DNAPositronium(const G4String& processName ="G4DNAPositronium");
  ~G4DNAPositronium();

  G4bool IsApplicable(const G4ParticleDefinition&);
  virtual void InitialiseProcess(const G4ParticleDefinition*);
  virtual void PrintInfo();
  
private:
  G4bool       isInitialised;

};


#endif
