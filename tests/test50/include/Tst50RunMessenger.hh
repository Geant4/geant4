#ifndef Tst50RunMessenger_h 
#define  Tst50RunMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class Tst50RunAction;
class G4UIcmdWithAString;
class G4UIdirectory;

class Tst50RunMessenger: public G4UImessenger
{
 
public:
  Tst50RunMessenger(Tst50RunAction*);

  ~Tst50RunMessenger();
 
void SetNewValue(G4UIcommand*, G4String);

private:  
  
  Tst50RunAction* p_Run;
  
   G4UIdirectory*  RunDir;
  G4UIcmdWithAString* testCmd;

};
#endif
