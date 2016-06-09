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
// $Id: G4UIWo.hh,v 1.4 2001/07/11 10:01:18 gunter Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
#ifndef G4UIWo_h
#define G4UIWo_h 

#if defined(G4UI_BUILD_WO_SESSION) || defined(G4UI_USE_WO)

#include "G4UIsession.hh"

class G4WoMessenger;

class G4UIWo : public G4UIsession {
public:
  G4UIWo(int,char**);
  ~G4UIWo();
  G4UIsession* SessionStart();
  void Prompt(G4String);
  void SessionTerminate();
  void PauseSessionStart(G4String);
  G4int ReceiveG4cout(G4String);
  G4int ReceiveG4cerr(G4String);
private:
  void SecondaryLoop(G4String);
  G4WoMessenger* woMessenger;
};

#endif

#endif

