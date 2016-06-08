
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The run messenger is defined
// Class Description - end
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em6RunMessenger_h
#define Em6RunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6RunAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6RunMessenger: public G4UImessenger
{
public: // Without description

   Em6RunMessenger(Em6RunAction* );
  ~Em6RunMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

private:

   Em6RunAction*              runAction;
   
   G4UIdirectory*             plotDir;   

#ifndef G4NOHIST
   G4UIcmdWithAString*        sethistNameCmd;
#endif

   G4UIcmdWithAnInteger*      setnbinStepCmd; 
   G4UIcmdWithADouble*        setSteplowCmd; 
   G4UIcmdWithADouble*        setStephighCmd; 

   G4UIcmdWithAnInteger*      setnbinEnCmd; 
   G4UIcmdWithADoubleAndUnit* setEnlowCmd; 
   G4UIcmdWithADoubleAndUnit* setEnhighCmd; 

   G4UIcmdWithAnInteger*      setnbinGammaCmd; 
   G4UIcmdWithADoubleAndUnit* setElowGammaCmd; 
   G4UIcmdWithADoubleAndUnit* setEhighGammaCmd; 

   G4UIcmdWithAnInteger*      setnbinTtCmd; 
   G4UIcmdWithADoubleAndUnit* setTtlowCmd; 
   G4UIcmdWithADoubleAndUnit* setTthighCmd; 

   G4UIcmdWithAnInteger*      setnbinTbCmd; 
   G4UIcmdWithADoubleAndUnit* setTblowCmd; 
   G4UIcmdWithADoubleAndUnit* setTbhighCmd; 

   G4UIcmdWithAnInteger*      setnbinTsecCmd; 
   G4UIcmdWithADoubleAndUnit* setTseclowCmd; 
   G4UIcmdWithADoubleAndUnit* setTsechighCmd; 

   G4UIcmdWithAnInteger*      setnbinRCmd; 
   G4UIcmdWithADoubleAndUnit* setRlowCmd; 
   G4UIcmdWithADoubleAndUnit* setRhighCmd; 

   G4UIcmdWithAnInteger*      setnbinThCmd; 
   G4UIcmdWithADoubleAndUnit* setThlowCmd; 
   G4UIcmdWithADoubleAndUnit* setThhighCmd; 

   G4UIcmdWithAnInteger*      setnbinThbackCmd; 
   G4UIcmdWithADoubleAndUnit* setThlowbackCmd; 
   G4UIcmdWithADoubleAndUnit* setThhighbackCmd; 

   G4UIcmdWithAnInteger*      setnbinzvertexCmd; 
   G4UIcmdWithADoubleAndUnit* setzlowCmd; 
   G4UIcmdWithADoubleAndUnit* setzhighCmd; 
 
};

#endif

