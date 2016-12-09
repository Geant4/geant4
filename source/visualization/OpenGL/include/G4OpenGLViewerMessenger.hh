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
//
// $Id: G4OpenGLViewerMessenger.hh 99440 2016-09-22 08:34:04Z gcosmo $

#ifndef G4OPENGLVIEWERMESSENGER_HH
#define G4OPENGLVIEWERMESSENGER_HH

#include "G4UImessenger.hh"

#include "G4String.hh"

class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class G4OpenGLViewerMessenger: public G4UImessenger {
public:
  static G4OpenGLViewerMessenger* GetInstance();  // Singleton constructor.
  ~G4OpenGLViewerMessenger();
  void SetNewValue (G4UIcommand*, G4String);

private:
  G4OpenGLViewerMessenger();  // Private constructor.
  static G4OpenGLViewerMessenger* fpInstance;
  G4UIdirectory* fpDirectory;
  G4UIcmdWithoutParameter* fpCommandPrintEPS;
  G4UIcommand* fpCommandPrintSize;
  G4UIcommand* fpCommandPrintFilename;
  G4UIcommand* fpCommandExport;
  G4UIcommand* fpCommandExportFormat;
  G4UIdirectory* fpDirectorySet;
  G4UIcommand* fpCommandDisplayHeadTime;
  G4UIcmdWithAnInteger* fpCommandDisplayListLimit;
  G4UIcommand* fpCommandDisplayLightFront;
  G4UIcommand* fpCommandEndTime;
  G4UIcmdWithAnInteger* fpCommandEventsDrawInterval;
  G4UIcmdWithADouble* fpCommandFade;
  G4UIcommand* fpCommandFlushAt;
  G4UIcmdWithAString* fpCommandPrintMode;
  G4UIcommand* fpCommandStartTime;
  G4UIcmdWithABool* fpCommandTransparency;
};

#endif
