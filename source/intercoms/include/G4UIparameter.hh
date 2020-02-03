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
//
// 
// ---------------------------------------------------------------------

#ifndef G4UIparameter_h
#define G4UIparameter_h 1

#include "globals.hh"
#include "G4UItokenNum.hh"

// class description:
//
//  This class represents a parameter which will be taken by a G4UIcommand
// object. In case a command is defined by constructing G4UIcmdXXX class,
// it automatically creates necessary parameter objects, thus the user needs 
// not to create parameter object(s) by him/herself. In case the user wants 
// to create a command directly instansiated by G4UIcommand class, he/she
// must create parameter object(s) by him/herself.

class G4UIparameter 
{
  public: // with description
      G4UIparameter();
      G4UIparameter(char theType);
      G4UIparameter(const char * theName, char theType, G4bool theOmittable);
      // Constructors, where "theName" is the name of the parameter which will
      // be used by the range checking, "theType" is the type of the parameter
      // (currently "b" (boolean), "i" (integer), "d" (double), and "s" (string)
      // are supported), and "theOmittable" is a boolean flag to set whether
      // the user of the command can ommit the parameter or not. If "theOmittable"
      // is true, the default value must be given.
      ~G4UIparameter();
      // Destructor. When a command is destructed, the delete operator(s) for 
      // associating parameter(s) are AUTOMATICALLY invoked. Thus the user needs 
      // NOT to invoke this by him/herself.

  public:
      G4bool operator==(const G4UIparameter &right) const;
      G4bool operator!=(const G4UIparameter &right) const;

      G4int CheckNewValue(const char* newValue);
      void List();

  private:
      G4String parameterName;
      G4String parameterGuidance;
      G4String defaultValue;
      G4String parameterRange;
      G4String parameterCandidate;
      char parameterType;
      G4bool omittable;
      G4bool currentAsDefaultFlag;
      G4int widget;

  public: // with description
      inline void SetDefaultValue(const char * theDefaultValue)
      { defaultValue = theDefaultValue; }
      void SetDefaultValue(G4int theDefaultValue);
      void SetDefaultValue(G4double theDefaultValue);
      // These methods set the default value of the parameter.
      void SetDefaultUnit(const char * theDefaultUnit);
      // This method can be used for a string-type parameter that is
      // used to specify a unit. This method is valid only for a
      // string-type parameter. With this set-method, not only the
      // default unit but also candidate units that belong to the
      // same unit category (a.k.a. dimension) as the default unit.
  public:
      inline G4String GetDefaultValue() const
      { return defaultValue; }
      inline char GetParameterType() const
      { return parameterType; }

  public: // with description
      inline void SetParameterRange(const char * theRange)
      { parameterRange = theRange; }
      //  Defines the range the parameter can take.
      //  The variable name appear in the range expression must be same
      // as the name of the parameter.
      //  All the C++ syntax of relational operators are allowed for the
      // range expression.
  public:
      inline G4String GetParameterRange() const
      { return parameterRange; }
    
    // parameterName
      inline void SetParameterName(const char * theName)
      { parameterName = theName; }
      inline G4String GetParameterName() const
      { return parameterName; }
    
  public: // with description
      inline void SetParameterCandidates(const char * theString)
      { parameterCandidate = theString; }
      //  This method is meaningful if the type of the parameter is string.
      // The candidates listed in the argument must be separated by space(s).
  public:
      inline G4String GetParameterCandidates() const
      { return parameterCandidate; }
    
    // omittable
      inline void SetOmittable(G4bool om)
      { omittable = om; }
      inline G4bool IsOmittable() const
      { return omittable; }
    
    // currentAsDefaultFlag
      inline void SetCurrentAsDefault(G4bool val)
      { currentAsDefaultFlag = val; }
      inline G4bool GetCurrentAsDefault() const
      { return currentAsDefaultFlag; }
    
    // out of date methods
      inline void SetWidget(G4int theWidget)
      { widget = theWidget; }
      inline const G4String GetParameterGuidance() const
      { return parameterGuidance; }
      inline void SetGuidance(const char * theGuidance)
      { parameterGuidance = theGuidance; }

  protected:
    using yystype = G4UItokenNum::yystype;
    using tokenNum = G4UItokenNum::tokenNum;
  private:
    // --- the following is used by CheckNewValue() -------
    G4int TypeCheck(const char* newValue );
    G4int RangeCheck(const char* newValue );
    G4int CandidateCheck(const char* newValue );
    G4int IsInt(const char* str, short maxDigit);
    G4int IsDouble(const char* str);
    G4int ExpectExponent(const char* str);
    //  syntax nodes
    yystype Expression( void );
    yystype LogicalORExpression( void );
    yystype LogicalANDExpression( void );
    yystype EqualityExpression ( void );
    yystype RelationalExpression( void );
    yystype AdditiveExpression( void );
    yystype MultiplicativeExpression( void );
    yystype UnaryExpression( void );
    yystype PrimaryExpression( void );
    //  semantics routines
    G4int Eval2( yystype arg1, G4int op, yystype arg2 );
    G4int CompareInt( G4int arg1, G4int op, G4int arg2);
    G4int CompareDouble( double arg1, G4int op, double arg2);
    //  utility 
    tokenNum Yylex( void );     // returns next token
    G4int G4UIpGetc( void );     // read one char from rangeBuf
    G4int G4UIpUngetc( G4int c );  // put back  
    G4int Backslash( G4int c );
    G4int Follow( G4int expect, G4int ifyes, G4int ifno );
    G4String TokenToStr(G4int token);
    //void PrintToken(void);  // debug
    //  data
    G4String rangeBuf;
    G4int bp;                  // buffer pointer for rangeBuf
    tokenNum token;
    yystype yylval;
    yystype newVal;
    G4int paramERR;
   //------------ end of CheckNewValue() related member --------------

};

#endif

