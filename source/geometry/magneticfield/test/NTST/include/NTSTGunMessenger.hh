//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: NTSTGunMessenger.hh,v 1.3 2006-06-29 18:25:38 gunter Exp $
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


