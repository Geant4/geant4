// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticlePropertyMessenger.hh,v 1.1 1999-01-07 16:10:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4ParticlePropertyMessenger.hh
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between ParticleDefinition and UI.
//
//  /particle/property/   Paricle Table control commands.
//   Commands : 
//     dump * dump particle properties.
//     stable * Set stable flag.
//     lifetime * Set life time.
//     Verbose * Set Verbose level
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    13 Nov. 1997, H. Kurashige   : fix bugs
//    08 Jan. 1998, H. Kurashige   : new UIcommand
//
//---------------------------------------------------------------

#ifndef G4ParticlePropertyMessenger_h
#define G4ParticlePropertyMessenger_h 1

class G4ParticleTable;
class G4ParticleDefinition;
class G4DecayTableMessenger;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger; 

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticlePropertyMessenger: public G4UImessenger
{
  public:
    G4ParticlePropertyMessenger(G4ParticleTable* pTable = NULL);
    ~G4ParticlePropertyMessenger();

    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4ParticlePropertyMessenger(const G4ParticlePropertyMessenger&){};

  private:
    G4ParticleDefinition* SetCurrentParticle();
    G4ParticleTable* theParticleTable;
    G4ParticleDefinition* currentParticle;

    G4UIdirectory *             thisDirectory;
    G4UIcmdWithoutParameter *   dumpCmd;
    G4UIcmdWithABool *          stableCmd; 
    G4UIcmdWithAnInteger *      verboseCmd;
    G4UIcmdWithADoubleAndUnit * lifetimeCmd; 
 
    G4DecayTableMessenger* fDecayTableMessenger;
};

#endif






