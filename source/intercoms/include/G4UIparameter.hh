// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIparameter.hh,v 1.2 1999-10-29 06:06:45 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ---------------------------------------------------------------------

#ifndef G4UIparameter_h
#define G4UIparameter_h 1

#include "globals.hh"

#ifdef WIN32
#  include <Strstrea.h>
#else
#  include <strstream.h>
#endif

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
      int operator==(const G4UIparameter &right) const;
      int operator!=(const G4UIparameter &right) const;

      G4int CheckNewValue(G4String newValue);
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

  private:
    // --- the following is used by CheckNewValue() -------
    int TypeCheck(G4String newValue );
    int RangeCheck(G4String newValue );
    int CandidateCheck(G4String newValue );
    int IsInt(const char* str, short maxDigit);
    int IsDouble(const char* str);
    int ExpectExponent(const char* str);
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
    int Eval2( yystype arg1, int op, yystype arg2 );
    int CompareInt( int arg1, int op, int arg2);
    int CompareDouble( double arg1, int op, double arg2);
    //  utility 
    tokenNum Yylex( void );     // returns next token
    int G4UIpGetc( void );     // read one char from rangeBuf
    int G4UIpUngetc( int c );  // put back  
    int Backslash( int c );
    int Follow( int expect, int ifyes, int ifno );
    G4String TokenToStr(int token);
    //void PrintToken(void);  // debug
    //  data
    G4String rangeBuf;
    int bp;                  // buffer pointer for rangeBuf
    tokenNum token;
    yystype yylval;
    yystype newVal;
    int paramERR;
   //------------ end of CheckNewValue() related member --------------

};

#endif

