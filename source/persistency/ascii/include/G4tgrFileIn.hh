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
// $Id: G4tgrFileIn.hh 68052 2013-03-13 14:38:53Z gcosmo $
//
//
// class G4tgrFileIn
//
// Class description:
//
// Singleton for importing file descriptions.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrFileIn_HH
#define G4tgrFileIn_HH

#include "globals.hh"
#include <vector>

class G4tgrFileIn 
{
  public:  // with description

    G4tgrFileIn();
   ~G4tgrFileIn();

    static G4tgrFileIn& GetInstance( const G4String& name ); 
      // Get the only instance opening the file

    static G4tgrFileIn& GetInstanceOpened( const G4String& name ); 
      // Get the only instance when file should be already opened
  
    G4int GetWordsInLine( std::vector<G4String>& wl );
      // Read a line and transform it to a vector of words 

    void ErrorInLine();
      // Print out an error message indicating the line being read

    // Access data members

    G4int Nline() { return theLineNo[theCurrentFile]; }
  
    const G4String& GetName() { return theName; }
  
    void OpenNewFile( const char* filename );
    G4bool EndOfFile();
    void Close();
    void DumpException( const G4String& sent );

  private:

    G4tgrFileIn( const G4String& name ) : theCurrentFile(-1), theName(name) {}

  private:

    std::vector< std::ifstream* > theFiles;

    std::vector<G4int> theLineNo;
      // Number of line being read

    std::vector<G4String> theNames;

    G4int theCurrentFile;
      // Index of file being read in theFiles
  
    static G4ThreadLocal std::vector<G4tgrFileIn*> *theInstances;
      // Vector of class instances (each one identified by its name)
  
    G4String theName;
      // Name of file
};

#endif 
