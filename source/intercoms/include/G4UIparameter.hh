// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UIparameter.hh,v 1.1 1999-01-07 16:09:23 gunter Exp $
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

class G4UIparameter 
{
  public:
      G4UIparameter();
      G4UIparameter(char theType);
      G4UIparameter(const char * theName, char theType, G4bool theOmittable);
      ~G4UIparameter();
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

  public:
    // defaultValue
      inline void SetDefaultValue(const char * theDefaultValue)
      { defaultValue = theDefaultValue; }
      void SetDefaultValue(G4int theDefaultValue);
      void SetDefaultValue(G4double theDefaultValue);
      inline G4String GetDefaultValue() const
      { return defaultValue; }

    // parameterType
      inline char GetParameterType() const
      { return parameterType; }

    // parameterRange
      inline void SetParameterRange(const char * theRange)
      { parameterRange = theRange; }
      inline G4String GetParameterRange() const
      { return parameterRange; }
    
    // parameterName
      inline void SetParameterName(const char * theName)
      { parameterName = theName; }
      inline G4String GetParameterName() const
      { return parameterName; }
    
    // parameterCandidates
      inline void SetParameterCandidates(const char * theString)
      { parameterCandidate = theString; }
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

