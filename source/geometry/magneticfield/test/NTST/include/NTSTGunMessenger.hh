//
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
// $Id: NTSTGunMessenger.hh,v 1.2 2003-12-09 15:35:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef NTSTGunMessenger_h
#define NTSTGunMessenger_h 1

class NTSTGunGenerator;
class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;

#include "G4UImessenger.hh"
#include "globals.hh"

// class description:
//
//  This is a concrete class of G4UImessenger which handles commands for
//  NTSTGunGenerator.
//

class NTSTGunMessenger: public G4UImessenger
{
public:
  NTSTGunMessenger(NTSTGunGenerator * NTSTGun);
  ~NTSTGunMessenger();
  
public:
  void SetNewValue(G4UIcommand * command,G4String newValues);
  G4String GetCurrentValue(G4UIcommand * command);
  
private:
  NTSTGunGenerator * fNTSTGun;
  G4ParticleTable * particleTable;
  
private: //commands
  G4UIdirectory *             gunDirectory;
  G4UIcmdWithoutParameter *   listCmd;
  G4UIcmdWithAString *        particleCmd;
  G4UIcmdWithADoubleAndUnit * plowCmd;
  G4UIcmdWithADoubleAndUnit * phighCmd;
  G4UIcmdWithADoubleAndUnit * t0Cmd;
  G4UIcmdWith3Vector *        polCmd;
  G4UIcmdWithAnInteger *      numberCmd;
  G4UIcmdWith3VectorAndUnit * meanVertexCmd;
  G4UIcmdWith3VectorAndUnit * rmsVertexCmd;
  G4UIcmdWithADouble *        coslowCmd;
  G4UIcmdWithADouble *        coshighCmd;
  G4UIcmdWithADoubleAndUnit * philowCmd;
  G4UIcmdWithADoubleAndUnit * phihighCmd;
};

#endif


