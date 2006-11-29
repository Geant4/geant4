#ifndef G4tgrUtils_h
#define G4tgrUtils_h 

#include <iostream>
using namespace std;
#include "globals.hh"

/*---------------------------------------------------------------------------
ClassName:   G4tgrUtils
Author:      P. Arce
Changes:     12/07/00: creation  
---------------------------------------------------------------------------*/ 
// Description  
//----------------------------------------------- 
/*! Class with several utilities, specially for managing of G4Strings */

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include <vector>
#include "CLHEP/Evaluator/Evaluator.h"

enum WLSIZEtype{WLSIZE_EQ,WLSIZE_NE,WLSIZE_LE,WLSIZE_LT,WLSIZE_GE,WLSIZE_GT};

class G4tgrUtils
{
public:
  G4tgrUtils(){};
  ~G4tgrUtils(){};
  
  //! Checks that every character in a G4String is a number (also '.' or exponencial notation: 'E')
  static bool IsSeparator(char);
  static bool IsNumber( const G4String& str);

  //! dumps a three vector with a message
  static void Dump3v( const CLHEP::Hep3Vector& vec, const char* msg);
  //! dumps a rotation matrix with a message
  static void Dumprm( const CLHEP::HepRotation& rm, const char* msg);
  //! dumps a vector of G4Strings with a message to cout
  static void DumpVS( const vector<G4String>& wl , const char* msg);
  //! dumps a vector of G4Strings with a message to outs
  static void DumpVS( const vector<G4String>& wl , const char* msg, ostream& outs) ;
  static G4String ftoa( const G4double dou );
  static void CheckWLsize( const vector<G4String>& wl, uint noWords, WLSIZEtype st, const G4String& methodName );


  //! Return the str without leading ':'
  static G4String SubColon( const G4String& str );

  //! Return the str without leading and trailing '"' and ' '
  static G4String SubQuotes( const G4String& str );

  //! Convert a G4String to an float, checking that it is really a number
  static double GetFloat( const G4String& str, double unitval = 1. );
  //! Convert a G4String to an integer, checking that it is really an integer
  static int GetInt( const G4String& str );
  //! Convert a bool to an integer, checking that it is really a bool
   static bool GetBool( const G4String& str );

private:
  static HepTool::Evaluator* theCLHEPevaluator;
};


#endif 


