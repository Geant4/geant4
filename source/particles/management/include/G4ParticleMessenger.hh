// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ParticleMessenger.hh,v 1.3 1999-10-28 23:24:12 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//---------------------------------------------------------------
//
//  G4ParticleMessenger.hh
//
//  Class Description:
//    This is a messenger class to interface to exchange information
//    between Particle related classes and UI.
//
// ------------------------------------------------------
//  the List of Directory and Commands
// ------------------------------------------------------
// G4ParticleMessenger
//  /particle/   Paricle control commands.
//   Commands : 
//    select * Select particle 
//    List * List name of particles.
//
// G4ParticlePropertyMessenger
//  /particle/property/   Paricle Table control commands.
//   Commands : 
//     dump * dump particle properties.
//     stable * Set stable flag.
//     lifetime * Set life time.
//     Verbose * Set Verbose level
//
// G4DecayTableMessenger
//  /particle/property/decay/   Decay Table control commands.
//   Commands : 
//     select * Enter index of decay mode.
//     dump * Dump decay mode information.
//     br * Set branching ratio. [0< BR <1.0]
//
//
//  History:
//    13 June 1997, H. Kurashige   : The 1st version created.
//    10 Nov 1997,  H.Kurashige    : add /particle/property/Verbose 
//    08 Jan. 1998, H. Kurashige   : new UIcommand
//    08 June 1998, H. Kurashige   : remove fProcessManagerMessenger
//    25 Nov. 1998, H. Kurashige   : add /particle/find
//---------------------------------------------------------------

#ifndef G4ParticleMessenger_h
#define G4ParticleMessenger_h 1

class G4ParticleDefinition;
class G4ParticleTable;
class G4ParticlePropertyMessenger;


class G4UIdirectory;
class G4UIcmdWithAString; 
class G4UIcmdWithAnInteger; 

#include "G4UImessenger.hh"
#include "globals.hh"


class G4ParticleMessenger: public G4UImessenger
{
  public: 
    G4ParticleMessenger(G4ParticleTable* pTable = 0);
    virtual ~G4ParticleMessenger();

  public: // With Description
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
  //  !!!  can not use "copy constructor" !!!!
    G4ParticleMessenger(const G4ParticleMessenger&){};

  private:
    G4UIdirectory *             thisDirectory;
    G4UIcmdWithAString *        listCmd;
    G4UIcmdWithAString *        selectCmd;
    G4UIcmdWithAnInteger *      findCmd;

    G4ParticleTable* theParticleTable;
    G4ParticleDefinition* currentParticle;

    G4ParticlePropertyMessenger* fParticlePropertyMessenger;
};

#endif









