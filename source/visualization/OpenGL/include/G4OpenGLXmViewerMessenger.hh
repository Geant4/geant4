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
//
// $Id: G4OpenGLXmViewerMessenger.hh,v 1.3 2006-01-11 18:42:59 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLXMVIEWERMESSENGER_HH
#define G4OPENGLXMVIEWERMESSENGER_HH

#include "G4UImessenger.hh"

#include "G4String.hh"

class G4OpenGLXmViewer;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

class G4OpenGLXmViewerMessenger: public G4UImessenger {
public:
  static G4OpenGLXmViewerMessenger* GetInstance();  // Singleton constructor.
  ~G4OpenGLXmViewerMessenger();
  void SetNewValue (G4UIcommand*, G4String);

private:
  G4OpenGLXmViewerMessenger();  // Private constructor.
  static G4OpenGLXmViewerMessenger* fpInstance;
  G4UIdirectory* fpDirectory;
  G4UIdirectory* fpDirectorySet;
  G4UIcmdWithADoubleAndUnit* fpCommandSetDollyHigh;
  G4UIcmdWithADoubleAndUnit* fpCommandSetDollyLow;
  G4UIcmdWithADoubleAndUnit* fpCommandSetPanHigh;
  G4UIcmdWithADoubleAndUnit* fpCommandSetRotationHigh;
  G4UIcmdWithADouble* fpCommandSetZoomHigh;
  G4UIcmdWithADouble* fpCommandSetZoomLow;
};

#endif

#endif

