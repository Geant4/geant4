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
// $Id: G4UIbatch.hh,v 1.6 2002-06-07 17:37:44 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $id$

#ifndef G4UIbatch_h
#define G4UIbatch_h 1

#include "globals.hh"
#include "G4UIsession.hh"
#include "g4std/fstream"

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
      G4std::ifstream macroFile;
      G4String macroFileName;
      G4bool openFailed;
};



#endif

