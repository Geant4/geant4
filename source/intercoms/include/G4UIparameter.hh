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
// G4UIparameter
//
// Class description:
//
// This class represents a parameter which will be taken by a G4UIcommand
// object. In case a command is defined by constructing G4UIcmdXXX class,
// it automatically creates necessary parameter objects, thus the user needs
// not to create parameter object(s). In case the user wants to create a
// command directly instantiated by G4UIcommand class, he/she must create
// a parameter object(s)

// Author: Makoto Asai, 1997
// --------------------------------------------------------------------
#ifndef G4UIparameter_hh
#define G4UIparameter_hh 1

#include "G4UItokenNum.hh"
#include "globals.hh"

class G4UIparameter
{
  public:
    // Default constructor
    G4UIparameter() = default;

    // Constructors, where "theName" is the name of the parameter which will
    // be used by the range checking, "theType" is the type of the parameter
    // (currently "b" (Boolean), "i" (integer), "l" (long int), "d" (double)
    // and "s" (string) are supported).
    // "theOmittable" is a Boolean flag to set whether
    // the user of the command can omit the parameter or not.
    // If "theOmittable" is true, the default value must be given
    G4UIparameter(char theType);
    G4UIparameter(const char* theName, char theType, G4bool theOmittable);

    // Destructor. When a command is destructed, the delete operator(s) of the
    // associated parameter(s) are AUTOMATICALLY invoked
    ~G4UIparameter();

    G4int CheckNewValue(const char* newValue);
    void List();

    // These methods set the default value of the parameter
    inline void SetDefaultValue(const char* theDefaultValue) { defaultValue = theDefaultValue; }
    void SetDefaultValue(G4int theDefaultValue);
    void SetDefaultValue(G4long theDefaultValue);
    void SetDefaultValue(G4double theDefaultValue);

    // This method can be used for a string-type parameter that is
    // used to specify a unit. This method is valid only for a
    // string-type parameter
    void SetDefaultUnit(const char* theDefaultUnit);

    inline const G4String& GetDefaultValue() const { return defaultValue; }
    inline char GetParameterType() const { return parameterType; }

    // Defines the range the parameter can take.
    // The variable name appearing in the range expression must be the
    // same as the name of the parameter.
    // All the C++ syntax of relational operators are allowed for the
    // range expression
    inline void SetParameterRange(const char* theRange) { rangeExpression = theRange; }

    inline const G4String& GetParameterRange() const { return rangeExpression; }

    inline void SetParameterName(const char* pName) { parameterName = pName; }
    inline const G4String& GetParameterName() const { return parameterName; }

    // This method is meaningful if the type of the parameter is string.
    // The candidates listed in the argument must be separated by space(s)
    inline void SetParameterCandidates(const char* theString) { parameterCandidate = theString; }

    inline const G4String& GetParameterCandidates() const { return parameterCandidate; }

    inline void SetOmittable(G4bool om) { omittable = om; }
    inline G4bool IsOmittable() const { return omittable; }

    inline void SetCurrentAsDefault(G4bool val) { currentAsDefaultFlag = val; }
    inline G4bool GetCurrentAsDefault() const { return currentAsDefaultFlag; }

    // Obsolete methods
    //
    inline const G4String& GetParameterGuidance() const { return parameterGuidance; }
    inline void SetGuidance(const char* theGuidance) { parameterGuidance = theGuidance; }

  protected:
    using yystype = G4UItokenNum::yystype;
    using tokenNum = G4UItokenNum::tokenNum;

  private:
    // --- the following is used by CheckNewValue() -------
    G4bool TypeCheck(const char* newValue);
    G4bool RangeCheck(const char* newValue);
    G4bool CandidateCheck(const char* newValue);
    //  syntax nodes
    yystype Expression();
    yystype LogicalORExpression();
    yystype LogicalANDExpression();
    yystype EqualityExpression();
    yystype RelationalExpression();
    yystype AdditiveExpression();
    yystype MultiplicativeExpression();
    yystype UnaryExpression();
    yystype PrimaryExpression();
    //  semantics routines
    G4int Eval2(const yystype& arg1, G4int op, const yystype& arg2);
    //  utility
    tokenNum Yylex();  // returns next token
    G4int G4UIpGetc();  // read one char from rangeBuf
    G4int G4UIpUngetc(G4int c);  // put back
    G4int Backslash(G4int c);
    G4int Follow(G4int expect, G4int ifyes, G4int ifno);

    // data -----------------------------------------------------------

    G4String parameterName;
    G4String parameterGuidance;
    G4String defaultValue;
    G4String rangeExpression;
    G4String parameterCandidate;
    char parameterType = '\0';
    G4bool omittable = false;
    G4bool currentAsDefaultFlag = false;

    //------------ CheckNewValue() related data members ---------------
    G4int bp = 0;  // current index in rangeExpression
    tokenNum token = G4UItokenNum::NONE;
    yystype yylval;
    yystype newVal;
    G4int paramERR = 0;
};

#endif
