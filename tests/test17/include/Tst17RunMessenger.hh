
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17RunMessenger.hh,v 1.1 1999-11-30 18:01:52 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

#ifndef Tst17RunMessenger_h
#define Tst17RunMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class Tst17RunAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

class Tst17RunMessenger: public G4UImessenger
{
  public:

   Tst17RunMessenger(Tst17RunAction* );
   ~Tst17RunMessenger();

   void SetNewValue(G4UIcommand* ,G4String );

  private:

   Tst17RunAction* runAction;
   
   G4UIdirectory* plotDir;   

   G4UIcmdWithAString* sethistNameCmd;

   G4UIcmdWithAnInteger* setnbinEnCmd; 
   G4UIcmdWithADoubleAndUnit* setEnlowCmd; 
   G4UIcmdWithADoubleAndUnit* setEnhighCmd; 

   G4UIcmdWithAnInteger* setnbinGammaCmd; 
   G4UIcmdWithADoubleAndUnit* setElowGammaCmd; 
   G4UIcmdWithADoubleAndUnit* setEhighGammaCmd; 

   G4UIcmdWithAnInteger* setnbinTtCmd; 
   G4UIcmdWithADoubleAndUnit* setTtlowCmd; 
   G4UIcmdWithADoubleAndUnit* setTthighCmd; 

   G4UIcmdWithAnInteger* setnbinTsecCmd; 
   G4UIcmdWithADoubleAndUnit* setTseclowCmd; 
   G4UIcmdWithADoubleAndUnit* setTsechighCmd; 
};

#endif

