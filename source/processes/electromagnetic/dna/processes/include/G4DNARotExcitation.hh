#ifndef G4DNARotExcitation_h
#define G4DNARotExcitation_h 1

#include "G4VEmProcess.hh"

class G4DNARotExcitation : public G4VEmProcess
{ 
public:
  G4DNARotExcitation(const G4String& processName ="G4DNARotExcitation");
  ~G4DNARotExcitation();

  G4bool IsApplicable(const G4ParticleDefinition&);
  virtual void InitialiseProcess(const G4ParticleDefinition*);
  virtual void PrintInfo();
  
private:
  G4bool       isInitialised;

};


#endif
