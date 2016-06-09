// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4ModelCommandsDrawByParticleID.hh,v 1.2 2005/11/23 20:24:15 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
// 
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description
// Trajectory model commands.
// Class Description - End:

#ifndef G4MODELCOMMANDDRAWBYPARTICLEID_HH
#define G4MODELCOMMANDDRAWBYPARTICLEID_HH

#include "G4String.hh"
#include "G4VModelCommand.hh"

class G4TrajectoryDrawByParticleID;
class G4UIcmdWithAString;
class G4UIcommand;

// Command to set positive/negative/neutral trajectory colouring through a string
class G4ModelCommandDrawByParticleIDSet : public G4VModelCommand<G4TrajectoryDrawByParticleID> {

public: // With description

  G4ModelCommandDrawByParticleIDSet(G4TrajectoryDrawByParticleID* model, const G4String& placement);

  virtual ~G4ModelCommandDrawByParticleIDSet();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithAString* fpCommand;

};

// Command to set default positive/negative/neutral trajectory colouring through G4Colour components
class G4ModelCommandDrawByParticleIDSetRGBA : public G4VModelCommand<G4TrajectoryDrawByParticleID> {

public:

  G4ModelCommandDrawByParticleIDSetRGBA(G4TrajectoryDrawByParticleID* model, const G4String& placement) ;
  virtual ~G4ModelCommandDrawByParticleIDSetRGBA();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithAString* fpCommand;

};

// Command to set default positive/negative/neutral trajectory colouring through a string
class G4ModelCommandDrawByParticleIDSetDefault : public G4VModelCommand<G4TrajectoryDrawByParticleID> {

public: // With description

  G4ModelCommandDrawByParticleIDSetDefault(G4TrajectoryDrawByParticleID* model, const G4String& placement);

  virtual ~G4ModelCommandDrawByParticleIDSetDefault();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithAString* fpCommand;

};

// Command to set default positive/negative/neutral trajectory colouring through G4Colour components
class G4ModelCommandDrawByParticleIDSetDefaultRGBA : public G4VModelCommand<G4TrajectoryDrawByParticleID> {

public:

  G4ModelCommandDrawByParticleIDSetDefaultRGBA(G4TrajectoryDrawByParticleID* model, const G4String& placement) ;
  virtual ~G4ModelCommandDrawByParticleIDSetDefaultRGBA();

  void SetNewValue(G4UIcommand* command, G4String newValue);

private:

  G4UIcmdWithAString* fpCommand;

};

#endif
