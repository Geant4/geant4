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
// $Id: G4tgrUtils.hh 66872 2013-01-15 01:25:57Z japost $
//
//
// class G4tgrUtils
//
// Class description:
//
// Utility class for managing of strings.

// History:
// - Created.                                 P.Arce, CIEMAT (November 2007)
// -------------------------------------------------------------------------

#ifndef G4tgrUtils_h
#define G4tgrUtils_h 

#include "globals.hh"

#include <iostream>
#include <vector>

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4tgrEvaluator.hh"

enum WLSIZEtype {WLSIZE_EQ,WLSIZE_NE,WLSIZE_LE,WLSIZE_LT,WLSIZE_GE,WLSIZE_GT};

class G4tgrUtils
{
  public:  // with description

    G4tgrUtils();
   ~G4tgrUtils();
  
    static G4bool IsSeparator(char);
    static G4bool IsNumber( const G4String& str);
    static G4bool IsInteger( const G4double val,
                             const G4double precision = 1.e-6 );
    static G4bool IsFunction( const G4String& word ); 
      // Checks that every character in a G4String is a number
      // (also '.' or exponencial notation: 'E')
    static G4bool WordIsUnit( const G4String& word );

    static void Dump3v( const G4ThreeVector& vec, const char* msg);
      // Dumps a three-vector with a message
    static void Dumprm( const G4RotationMatrix& rm, const char* msg);
      // Dumps a rotation matrix with a message
    static void DumpVS( const std::vector<G4String>& wl , const char* msg);
      // Dumps a vector of G4Strings with a message to cout
    static void DumpVS( const std::vector<G4String>& wl , const char* msg,
                              std::ostream& outs) ;
      // Dumps a vector of G4Strings with a message to outs
    static void CheckWLsize( const std::vector<G4String>& wl,
                                   unsigned int nWCheck, WLSIZEtype st,
                             const G4String& methodName );
    static G4bool CheckListSize( unsigned int nWreal,
                                 unsigned int nWcheck, WLSIZEtype st,
                                              G4String& outstr );

    static G4String SubColon( const G4String& str );
      // Return the str without leading ':'

    static G4String GetString( const G4String& str );
      // Return the str without leading and trailing '"' and ' '

    static G4double GetDouble( const G4String& str, G4double unitval = 1. );
      // Convert a G4String to a double, checking that it is really a number
    static G4int GetInt( const G4String& str );
      // Convert a G4String to an integer, checking that it is really an int
    static G4bool GetBool( const G4String& str );
      // Convert a bool to an integer, checking that it is really a bool

    static G4RotationMatrix GetRotationFromDirection( G4ThreeVector dir );

    static G4bool AreWordsEquivalent( const G4String& word1,
                                      const G4String& word2 );
      // Looks if word1 and word2 are equivalent, considering that
      // word1 may have '*', meaning 'any character'

  private:

    static G4ThreadLocal G4tgrEvaluator* theEvaluator;
};

#endif 
