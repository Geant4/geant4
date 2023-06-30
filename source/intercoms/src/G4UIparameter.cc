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
// Author: Makoto Asai, 1997
// --------------------------------------------------------------------

#include "G4UIparameter.hh"

#include "G4Tokenizer.hh"
#include "G4UIcommand.hh"
#include "G4UIcommandStatus.hh"
#include "G4UIparsing.hh"
#include "G4ios.hh"

#include <cctype>  // for CheckNewValue()
#include <sstream>

using namespace G4UItokenNum;

// --------------------------------------------------------------------
G4UIparameter::G4UIparameter(char theType)
{
  parameterType = theType;
}

// --------------------------------------------------------------------
G4UIparameter::G4UIparameter(const char* theName, char theType, G4bool theOmittable)
{
  parameterName = theName;
  parameterType = theType;
  omittable = theOmittable;
}

// --------------------------------------------------------------------
G4UIparameter::~G4UIparameter() = default;

// --------------------------------------------------------------------
void G4UIparameter::List()
{
  G4cout << G4endl << "Parameter : " << parameterName << G4endl;
  if (!parameterGuidance.empty()) {
    G4cout << parameterGuidance << G4endl;
  }
  G4cout << " Parameter type  : " << parameterType << G4endl;
  if (omittable) {
    G4cout << " Omittable       : True" << G4endl;
  }
  else {
    G4cout << " Omittable       : False" << G4endl;
  }
  if (currentAsDefaultFlag) {
    G4cout << " Default value   : taken from the current value" << G4endl;
  }
  else if (!defaultValue.empty()) {
    G4cout << " Default value   : " << defaultValue << G4endl;
  }
  if (!rangeExpression.empty()) {
    G4cout << " Parameter range : " << rangeExpression << G4endl;
  }
  if (!parameterCandidate.empty()) {
    G4cout << " Candidates      : " << parameterCandidate << G4endl;
  }
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultValue(G4int theDefaultValue)
{
  defaultValue = G4UIparsing::TtoS(theDefaultValue);
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultValue(G4long theDefaultValue)
{
  defaultValue = G4UIparsing::TtoS(theDefaultValue);
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultValue(G4double theDefaultValue)
{
  defaultValue = G4UIparsing::TtoS(theDefaultValue);
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultUnit(const char* theDefaultUnit)
{
  char type = (char)std::toupper(parameterType);
  if (type != 'S') {
    G4ExceptionDescription ed;
    ed << "This method can be used only for a string-type parameter that is "
          "used to specify a unit.\n"
       << "This parameter <" << parameterName << "> is defined as ";
    switch (type) {
      case 'D':
        ed << "double.";
        break;
      case 'I':
        ed << "integer.";
        break;
      case 'L':
        ed << "long int.";
        break;
      case 'B':
        ed << "bool.";
        break;
      default:
        ed << "undefined.";
    }
    G4Exception("G4UIparameter::SetDefaultUnit", "INTERCOM2010", FatalException, ed);
  }
  SetDefaultValue(theDefaultUnit);
  SetParameterCandidates(G4UIcommand::UnitsList(G4UIcommand::CategoryOf(theDefaultUnit)));
}

// ---------- CheckNewValue() related routines ------------------------

G4int G4UIparameter::CheckNewValue(const char* newValue)
{
  if (!TypeCheck(newValue)) {
    return fParameterUnreadable;
  }
  if (!RangeCheck(newValue)) {
    return fParameterOutOfRange;
  }
  if (!CandidateCheck(newValue)) {
    return fParameterOutOfCandidates;
  }
  return 0;  // succeeded
}

// --------------------------------------------------------------------
G4bool G4UIparameter::CandidateCheck(const char* newValue)
{
  if (parameterCandidate.empty()) {
    return true;
  }

  G4Tokenizer candidateTokenizer(parameterCandidate);
  G4String aToken;
  while (!(aToken = candidateTokenizer()).empty()) {
    if (aToken == newValue) {
      return true;
    }
  }
  G4cerr << "parameter value (" << newValue << ") is not listed in the candidate List." << G4endl;
  G4cerr << "  Candidates are:";
  G4Tokenizer candidateListTokenizer(parameterCandidate);
  while (!(aToken = candidateListTokenizer()).empty()) {
    G4cerr << ' ' << aToken;
  }
  G4cerr << G4endl;

  return false;
}

// --------------------------------------------------------------------
G4bool G4UIparameter::TypeCheck(const char* newValue)
{
  G4String newValueString(newValue);
  char type = (char)std::toupper(parameterType);
  switch (type) {
    case 'D':
      if (!G4UIparsing::IsDouble(newValueString.data())) {
        G4cerr << newValue << ": double value expected." << G4endl;
        return false;
      }
      break;
    case 'I':
      if (!G4UIparsing::IsInt(newValueString.data(), 10)) {
        G4cerr << newValue << ": integer expected." << G4endl;
        return false;
      }
      break;
    case 'L':
      if (!G4UIparsing::IsInt(newValueString.data(), 20)) {
        G4cerr << newValue << ": long int expected." << G4endl;
        return false;
      }
      break;
    case 'S':
      break;
    case 'B':
      G4StrUtil::to_upper(newValueString);
      if (newValueString == "Y" || newValueString == "N" || newValueString == "YES"
          || newValueString == "NO" || newValueString == "1" || newValueString == "0"
          || newValueString == "T" || newValueString == "F" || newValueString == "TRUE"
          || newValueString == "FALSE")
      {
        return true;
      }

      G4cerr << newValue << ": bool expected." << G4endl;
      return false;

    default:;
  }
  return true;
}

// --------------------------------------------------------------------
G4bool G4UIparameter::RangeCheck(const char* newValue)
{
  if (rangeExpression.empty()) {
    return true;
  }

  yystype result;
  bp = 0;  // reset buffer pointer for G4UIpGetc()
  std::istringstream is(newValue);
  char type = (char)std::toupper(parameterType);
  switch (type) {
    case 'D':
      is >> newVal.D;
      break;
    case 'I':
      is >> newVal.I;
      break;
    case 'L':
      is >> newVal.L;
      break;
    default:;
  }
  // PrintToken();          // Print tokens (consumes all tokens)
  token = Yylex();
  result = Expression();
  if (paramERR == 1) {
    return false;
  }
  if (result.type != CONSTINT) {
    G4cerr << "Illegal Expression in parameter range." << G4endl;
    return false;
  }
  if (result.I != 0) {
    return true;
  }
  G4cerr << "parameter out of range: " << rangeExpression << G4endl;
  return false;
}

// ------------------ syntax node functions  ------------------

yystype G4UIparameter::Expression()
{
  return LogicalORExpression();
}

// --------------------------------------------------------------------
yystype G4UIparameter::LogicalORExpression()
{
  yystype result;
  yystype p;
  p = LogicalANDExpression();
  if (token != LOGICALOR) {
    return p;
  }
  if (p.type == CONSTSTRING || p.type == IDENTIFIER) {
    G4cerr << "Parameter range: illegal type at '||'" << G4endl;
    paramERR = 1;
  }
  result.I = p.I;
  while (token == LOGICALOR) {
    token = Yylex();
    p = LogicalANDExpression();
    if (p.type == CONSTSTRING || p.type == IDENTIFIER) {
      G4cerr << "Parameter range: illegal type at '||'" << G4endl;
      paramERR = 1;
    }
    switch (p.type) {
      case CONSTINT:
        result.I += p.I;
        result.type = CONSTINT;
        break;
      case CONSTLONG:
        result.I += static_cast<int>(p.L != 0L);
        result.type = CONSTINT;
        break;
      case CONSTDOUBLE:
        result.I += static_cast<int>(p.D != 0.0);
        result.type = CONSTINT;
        break;
      default:
        G4cerr << "Parameter range: unknown type" << G4endl;
        paramERR = 1;
    }
  }
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::LogicalANDExpression()
{
  yystype result;
  yystype p;
  p = EqualityExpression();
  if (token != LOGICALAND) {
    return p;
  }
  if (p.type == CONSTSTRING || p.type == IDENTIFIER) {
    G4cerr << "Parameter range: illegal type at '&&'" << G4endl;
    paramERR = 1;
  }
  result.I = p.I;
  while (token == LOGICALAND) {
    token = Yylex();
    p = EqualityExpression();
    if (p.type == CONSTSTRING || p.type == IDENTIFIER) {
      G4cerr << "Parameter range: illegal type at '&&'" << G4endl;
      paramERR = 1;
    }
    switch (p.type) {
      case CONSTINT:
        result.I *= p.I;
        result.type = CONSTINT;
        break;
      case CONSTLONG:
        result.I *= static_cast<int>(p.L != 0L);
        result.type = CONSTINT;
        break;
      case CONSTDOUBLE:
        result.I *= static_cast<int>(p.D != 0.0);
        result.type = CONSTINT;
        break;
      default:
        G4cerr << "Parameter range: unknown type." << G4endl;
        paramERR = 1;
    }
  }
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::EqualityExpression()
{
  yystype arg1, arg2;
  G4int operat;
  yystype result;
  result = RelationalExpression();
  if (token == EQ || token == NE) {
    operat = token;
    token = Yylex();
    arg1 = result;
    arg2 = RelationalExpression();
    result.I = Eval2(arg1, operat, arg2);  // semantic action
    result.type = CONSTINT;
  }
  else {
    if (result.type != CONSTINT && result.type != CONSTDOUBLE) {
      G4cerr << "Parameter range: error at EqualityExpression" << G4endl;
      paramERR = 1;
    }
  }
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::RelationalExpression()
{
  yystype arg1, arg2;
  G4int operat;
  yystype result;

  arg1 = AdditiveExpression();
  if (token == GT || token == GE || token == LT || token == LE) {
    operat = token;
    token = Yylex();
    arg2 = AdditiveExpression();
    result.I = Eval2(arg1, operat, arg2);  // semantic action
    result.type = CONSTINT;
  }
  else {
    result = arg1;
  }
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::AdditiveExpression()
{
  yystype result = MultiplicativeExpression();
  if (token != '+' && token != '-') {
    return result;
  }
  G4cerr << "Parameter range: operator " << (char)token << " is not supported." << G4endl;
  paramERR = 1;
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::MultiplicativeExpression()
{
  yystype result = UnaryExpression();
  if (token != '*' && token != '/' && token != '%') {
    return result;
  }
  G4cerr << "Parameter range: operator " << (char)token << " is not supported." << G4endl;
  paramERR = 1;
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::UnaryExpression()
{
  yystype result;
  yystype p;
  switch (token) {
    case '-':
      token = Yylex();
      p = UnaryExpression();
      if (p.type == CONSTINT) {
        result.I = -p.I;
        result.type = CONSTINT;
      }
      if (p.type == CONSTLONG) {
        result.L = -p.L;
        result.type = CONSTLONG;
      }
      if (p.type == CONSTDOUBLE) {
        result.D = -p.D;
        result.type = CONSTDOUBLE;
      }
      break;
    case '+':
      token = Yylex();
      result = UnaryExpression();
      break;
    case '!':
      token = Yylex();
      G4cerr << "Parameter range error: "
             << "operator '!' is not supported (sorry)." << G4endl;
      paramERR = 1;
      result = UnaryExpression();
      break;
    default:
      result = PrimaryExpression();
  }
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::PrimaryExpression()
{
  yystype result;
  switch (token) {
    case IDENTIFIER:
      result.S = yylval.S;
      result.type = token;
      token = Yylex();
      break;
    case CONSTINT:
      result.I = yylval.I;
      result.type = token;
      token = Yylex();
      break;
    case CONSTLONG:
      result.L = yylval.L;
      result.type = token;
      token = Yylex();
      break;
    case CONSTDOUBLE:
      result.D = yylval.D;
      result.type = token;
      token = Yylex();
      break;
    case '(':
      token = Yylex();
      result = Expression();
      if (token != ')') {
        G4cerr << " ')' expected" << G4endl;
        paramERR = 1;
      }
      token = Yylex();
      break;
    default:
      return result;
  }
  return result;  // never executed
}

//---------------- semantic routines ---------------------------------

G4int G4UIparameter::Eval2(const yystype& arg1, G4int op, const yystype& arg2)
{
  if ((arg1.type != IDENTIFIER) && (arg2.type != IDENTIFIER)) {
    G4cerr << parameterName << ": meaningless comparison " << G4int(arg1.type) << " "
           << G4int(arg2.type) << G4endl;
    paramERR = 1;
  }
  char type = (char)std::toupper(parameterType);
  if (arg1.type == IDENTIFIER) {
    switch (type) {
      case 'I':
        if (arg2.type == CONSTINT) {
          return G4UIparsing::CompareInt(newVal.I, op, arg2.I, paramERR);
        }
        else {
          G4cerr << "integer operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'L':
        if (arg2.type == CONSTLONG) {
          return G4UIparsing::CompareLong(newVal.L, op, arg2.L, paramERR);
        }
        else {
          G4cerr << "long int operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'D':
        if (arg2.type == CONSTDOUBLE) {
          return G4UIparsing::CompareDouble(newVal.D, op, arg2.D, paramERR);
        }
        else if (arg2.type == CONSTINT) {  // integral promotion
          return G4UIparsing::CompareDouble(newVal.D, op, arg2.I, paramERR);
        }
        else if (arg2.type == CONSTLONG) {
          return G4UIparsing::CompareDouble(newVal.D, op, arg2.L, paramERR);
        }
        break;
      default:;
    }
  }
  if (arg2.type == IDENTIFIER) {
    switch (type) {
      case 'I':
        if (arg1.type == CONSTINT) {
          return G4UIparsing::CompareInt(arg1.I, op, newVal.I, paramERR);
        }
        else {
          G4cerr << "integer operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'L':
        if (arg1.type == CONSTLONG) {
          return G4UIparsing::CompareLong(arg1.L, op, newVal.L, paramERR);
        }
        else {
          G4cerr << "long int operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'D':
        if (arg1.type == CONSTDOUBLE) {
          return G4UIparsing::CompareDouble(arg1.D, op, newVal.D, paramERR);
        }
        else if (arg1.type == CONSTINT) {  // integral promotion
          return G4UIparsing::CompareDouble(arg1.I, op, newVal.D, paramERR);
        }
        else if (arg1.type == CONSTLONG) {  // integral promotion
          return G4UIparsing::CompareDouble(arg1.L, op, newVal.D, paramERR);
        }
        break;
      default:;
    }
  }
  G4cerr << "no param name is specified at the param range." << G4endl;
  return 0;
}

// --------------------- utility functions --------------------------

tokenNum G4UIparameter::Yylex()  // reads input and returns token number KR486
{  // (returns EOF)
  G4int c;
  G4String buf;

  while ((c = G4UIpGetc()) == ' ' || c == '\t' || c == '\n') {
    ;
  }
  if (c == EOF) {
    return (tokenNum)EOF;  // KR488
  }
  buf = "";
  if ((isdigit(c) != 0) || c == '.') {  // I or D
    do {
      buf += (unsigned char)c;
      c = G4UIpGetc();
    } while (c == '.' || (isdigit(c) != 0) || c == 'e' || c == 'E' || c == '+' || c == '-');
    G4UIpUngetc(c);
    const char* t = buf;
    std::istringstream is(t);
    if (G4UIparsing::IsInt(buf.data(), 20)) {
      is >> yylval.I;
      return CONSTINT;
    }
    if (G4UIparsing::IsDouble(buf.data())) {
      is >> yylval.D;
      return CONSTDOUBLE;
    }

    G4cerr << buf << ": numeric format error." << G4endl;
  }
  buf = "";
  if ((isalpha(c) != 0) || c == '_') {  // IDENTIFIER
    do {
      buf += (unsigned char)c;
    } while ((c = G4UIpGetc()) != EOF && ((isalnum(c) != 0) || c == '_'));
    G4UIpUngetc(c);
    if (buf == parameterName) {
      yylval.S = buf;
      return IDENTIFIER;
    }

    G4cerr << buf << " is not a parameter name." << G4endl;
    paramERR = 1;
  }
  switch (c) {
    case '>':
      return (tokenNum)Follow('=', GE, GT);
    case '<':
      return (tokenNum)Follow('=', LE, LT);
    case '=':
      return (tokenNum)Follow('=', EQ, '=');
    case '!':
      return (tokenNum)Follow('=', NE, '!');
    case '|':
      return (tokenNum)Follow('|', LOGICALOR, '|');
    case '&':
      return (tokenNum)Follow('&', LOGICALAND, '&');
    default:
      return (tokenNum)c;
  }
}

// --------------------------------------------------------------------
G4int G4UIparameter::Follow(G4int expect, G4int ifyes, G4int ifno)
{
  G4int c = G4UIpGetc();
  if (c == expect) {
    return ifyes;
  }
  G4UIpUngetc(c);
  return ifno;
}

//------------------ low level routines -----------------------------

G4int G4UIparameter::G4UIpGetc()
{  // emulation of getc()
  auto length = (G4int)rangeExpression.length();
  if (bp < length) {
    return rangeExpression[bp++];
  }

  return EOF;
}

// --------------------------------------------------------------------
G4int G4UIparameter::G4UIpUngetc(G4int c)
{  // emulation of ungetc()
  if (c < 0) {
    return -1;
  }
  if (bp > 0 && c == rangeExpression[bp - 1]) {
    --bp;
  }
  else {
    G4cerr << "G4UIpUngetc() failed." << G4endl;
    G4cerr << "bp=" << bp << " c=" << c << " pR(bp-1)=" << rangeExpression[bp - 1] << G4endl;
    paramERR = 1;
    return -1;
  }
  return 0;
}
