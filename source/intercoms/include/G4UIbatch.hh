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
// $Id: G4UIbatch.hh,v 1.9 2006/06/29 19:07:44 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// $id$

#ifndef G4UIbatch_h
#define G4UIbatch_h 1

#include "globals.hh"
#include "G4UIsession.hh"
#include <fstream>

class G4UImanager;

// class description
//  This is a concrete class of G4UIsession.
//  This class is constructed by G4UImanager when the user 
// applies "/control/execute macro_file" command.
//  If the user's application runs in pure batch mode with
// only one fixed-name macro file, he/she can construct this
// class object with giving the file name and use it as other
// ordinary G4UIsession concrete class, i.e. invoke SessionStart()
// from his/her main().

class G4UIbatch : public G4UIsession 
{
  public: // with description
      G4UIbatch(const char* fileName,G4UIsession* prevSession=NULL);
      //  Constructor. 
      //  "prevSession" must be NULL if this class is constructed
      // from main().
  public:
      ~G4UIbatch();

      G4UIsession * SessionStart();
      void PauseSessionStart(G4String Prompt);
  
  private:
      G4UImanager * UImanager;
      G4UIsession * previousSession;
      std::ifstream macroFile;
      G4String macroFileName;
      G4bool openFailed;

  private:
      static G4bool commandFailed;
};



#endif

