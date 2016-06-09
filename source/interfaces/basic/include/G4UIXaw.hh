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
// $Id: G4UIXaw.hh,v 1.6 2006-06-29 19:09:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4UIXaw_h
#define G4UIXaw_h 

#if defined(G4UI_BUILD_XAW_SESSION) || defined(G4UI_USE_XAW)

#include <X11/Intrinsic.h>

#include "G4VBasicShell.hh"

class G4UIsession;
class G4UImanager;

// Class description :
//
//  G4UIXaw : class to handle an interactive session.
// G4UIXaw is the "Athena" version of G4UIterminal.
// It has been provided at a time where no free Motif
// where available. Now that some exists (lesstif), 
// G4UIXm should be prefered. The look and feel of 
// G4UIXaw is hugly...
//
//  No command completion is available.
//
// Class description - end :

class G4UIXaw : public G4VBasicShell {
public:
  G4UIXaw(int,char**);
  ~G4UIXaw();
  G4UIsession* SessionStart();
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(G4String);
  Widget GetDialog();
private:
  void SecondaryLoop(G4String);
  G4bool GetHelpChoice(G4int&);
  void ExitHelp();
private:
  Widget shell,dialog;
  G4bool fHelp;
  G4int fHelpChoice;
  static void Callback(Widget,XtPointer,XtPointer);
};

#endif

#endif

