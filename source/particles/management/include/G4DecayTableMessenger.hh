// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4DecayTableMessenger.hh,v 1.1 1999-01-07 16:10:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4DecayTableMessenger.hh
//
//  Description:
//    This is a messenger class to interface to exchange information
//    between Decay Table/Decay Channel and UI.
//
// G4DecayTableMessenger
//  /particle/property/decay/   Decay Table control commands.
//   Commands : 
//     select * Enter index of decay mode.
//     dump * Dump decay mode information.
//     br * Set branching ratio. [0< BR <1.0]
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    13 Nov. 1997, H. Kurashige   : fix bugs
//    08 Jan. 1998, H. Kurashige   : new UIcommand
//---------------------------------------------------------------

#ifndef G4DecayTableMessenger_h
#define G4DecayTableMessenger_h 1

class G4ParticleTable;
class G4VDecayChannel;
class G4ParticleDefinition;
class G4DecayTable;

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAnInteger; 
class G4UIcmdWithADouble;

#include "G4UImessenger.hh"
#include "globals.hh"

class G4DecayTableMessenger: public G4UImessenger
{
  public:
    G4DecayTableMessenger(G4ParticleTable* pTable = NULL);
    ~G4DecayTableMessenger();

    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    G4DecayTableMessenger(const G4DecayTableMessenger&){};

  private:
    G4ParticleDefinition* SetCurrentParticle();
    G4ParticleTable* theParticleTable;
    G4ParticleDefinition* currentParticle;
    G4DecayTable*   currentDecayTable;
    G4int           idxCurrentChannel;
    G4VDecayChannel* currentChannel;

    G4UIdirectory *             thisDirectory;
    G4UIcmdWithoutParameter *   dumpCmd;
    G4UIcmdWithAnInteger *      selectCmd;
    G4UIcmdWithADouble   *      brCmd; 
};

#endif






