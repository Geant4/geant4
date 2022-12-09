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
#include "G4UIcommandStatus.hh"
#include "G4Tokenizer.hh"
#include "G4ios.hh"
#include "G4UIcommand.hh"

#include <sstream>
#include <ctype.h>  // for CheckNewValue()

using namespace G4UItokenNum;

// --------------------------------------------------------------------
G4UIparameter::G4UIparameter(char theType)
{
  parameterType        = theType;
}

// --------------------------------------------------------------------
G4UIparameter::G4UIparameter(const char* theName, char theType,
                             G4bool theOmittable)
{
  parameterName = theName;
  parameterType = theType;
  omittable     = theOmittable;
}

// --------------------------------------------------------------------
G4UIparameter::~G4UIparameter() = default;

// --------------------------------------------------------------------
G4bool G4UIparameter::operator==(const G4UIparameter& right) const
{
  return (this == &right);
}

// --------------------------------------------------------------------
G4bool G4UIparameter::operator!=(const G4UIparameter& right) const
{
  return (this != &right);
}

// --------------------------------------------------------------------
void G4UIparameter::List()
{
  G4cout << G4endl << "Parameter : " << parameterName << G4endl;
  if(!parameterGuidance.empty())
  {
    G4cout << parameterGuidance << G4endl;
  }
  G4cout << " Parameter type  : " << parameterType << G4endl;
  if(omittable)
  {
    G4cout << " Omittable       : True" << G4endl;
  }
  else
  {
    G4cout << " Omittable       : False" << G4endl;
  }
  if(currentAsDefaultFlag)
  {
    G4cout << " Default value   : taken from the current value" << G4endl;
  }
  else if(!defaultValue.empty())
  {
    G4cout << " Default value   : " << defaultValue << G4endl;
  }
  if(!parameterRange.empty())
  {
    G4cout << " Parameter range : " << parameterRange << G4endl;
  }
  if(!parameterCandidate.empty())
  {
    G4cout << " Candidates      : " << parameterCandidate << G4endl;
  }
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultValue(G4int theDefaultValue)
{
  std::ostringstream os;
  os << theDefaultValue;
  defaultValue = os.str();
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultValue(G4long theDefaultValue)
{
  std::ostringstream os;
  os << theDefaultValue;
  defaultValue = os.str();
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultValue(G4double theDefaultValue)
{
  std::ostringstream os;
  os << theDefaultValue;
  defaultValue = os.str();
}

// --------------------------------------------------------------------
void G4UIparameter::SetDefaultUnit(const char* theDefaultUnit)
{
  char type = (char)std::toupper(parameterType);
  if(type != 'S')
  {
    G4ExceptionDescription ed;
    ed << "This method can be used only for a string-type parameter that is "
          "used to specify a unit.\n"
       << "This parameter <" << parameterName << "> is defined as ";
    switch(type)
    {
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
    G4Exception("G4UIparameter::SetDefaultUnit", "INTERCOM2010", FatalException,
                ed);
  }
  SetDefaultValue(theDefaultUnit);
  SetParameterCandidates(
    G4UIcommand::UnitsList(G4UIcommand::CategoryOf(theDefaultUnit)));
}

// ---------- CheckNewValue() related routines ------------------------

//#include "checkNewValue_debug.icc"
//#define DEBUG 1

// --------------------------------------------------------------------
G4int G4UIparameter::CheckNewValue(const char* newValue)
{
  if(TypeCheck(newValue) == 0)
  {
    return fParameterUnreadable;
  }
  if(!parameterRange.empty())
  {
    if(RangeCheck(newValue) == 0)
    {
      return fParameterOutOfRange;
    }
  }
  if(!parameterCandidate.empty())
  {
    if(CandidateCheck(newValue) == 0)
    {
      return fParameterOutOfCandidates;
    }
  }
  return 0;  // succeeded
}

// --------------------------------------------------------------------
G4int G4UIparameter::CandidateCheck(const char* newValue)
{
  G4Tokenizer candidateTokenizer(parameterCandidate);
  G4String aToken;
  G4int iToken = 0;
  while(!(aToken = candidateTokenizer()).empty())
  {
    ++iToken;
    if(aToken == newValue)
    {
      return iToken;
    }
  }
  G4cerr << "parameter value (" << newValue
         << ") is not listed in the candidate List." << G4endl;
  G4cerr << "  Candidates are:";
  G4Tokenizer candidateListTokenizer(parameterCandidate);
  while(!(aToken = candidateListTokenizer()).empty())
  {
    G4cerr << ' ' << aToken;
  }
  G4cerr << G4endl;

  return 0;
}

// --------------------------------------------------------------------
G4int G4UIparameter::RangeCheck(const char* newValue)
{
  yystype result;
  bp = 0;  // reset buffer pointer for G4UIpGetc()
  std::istringstream is(newValue);
  char type = (char)std::toupper(parameterType);
  switch(type)
  {
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
  token  = Yylex();
  result = Expression();
  if(paramERR == 1)
  {
    return 0;
  }
  if(result.type != CONSTINT)
  {
    G4cerr << "Illegal Expression in parameter range." << G4endl;
    return 0;
  }
  if(result.I != 0)
  {
    return 1;
  }
  G4cerr << "parameter out of range: " << parameterRange << G4endl;
  return 0;
}

// --------------------------------------------------------------------
G4int G4UIparameter::TypeCheck(const char* newValue)
{
  G4String newValueString(newValue);
  char type = (char)std::toupper(parameterType);
  switch(type)
  {
    case 'D':
      if(IsDouble(newValueString.data()) == 0)
      {
        G4cerr << newValue << ": double value expected." << G4endl;
        return 0;
      }
      break;
    case 'I':
      if(IsInt(newValueString.data(), 10) == 0)
      {
        G4cerr << newValue << ": integer expected." << G4endl;
        return 0;
      }
      break;
    case 'L':
      if(IsInt(newValueString.data(), 20) == 0)
      {
        G4cerr << newValue << ": long int expected." << G4endl;
        return 0;
      }
      break;
    case 'S':
      break;
    case 'B':
      G4StrUtil::to_upper(newValueString);
      if(newValueString == "Y" || newValueString == "N" ||
         newValueString == "YES" || newValueString == "NO" ||
         newValueString == "1" || newValueString == "0" ||
         newValueString == "T" || newValueString == "F" ||
         newValueString == "TRUE" || newValueString == "FALSE")
      {
        return 1;
      }
      else
      {
        G4cerr << newValue << ": bool expected." << G4endl;
        return 0;
      }
    default:;
  }
  return 1;
}

// --------------------------------------------------------------------
G4int G4UIparameter::IsInt(const char* buf,
                           short maxDigits)  // do not allow any std::ws
{
  const char* p = buf;
  G4int length  = 0;
  if(*p == '+' || *p == '-')
  {
    ++p;
  }
  if(isdigit((G4int) (*p)) != 0)
  {
    while(isdigit((G4int) (*p)) != 0)
    {
      ++p;
      ++length;
    }
    if(*p == '\0')
    {
      if(length > maxDigits)
      {
        G4cerr << "digit length exceeds" << G4endl;
        return 0;
      }
      return 1;
    }
    else
    {
      // G4cerr <<"illegal character after int:"<<buf<<G4endl;
    }
  }
  else
  {
    // G4cerr <<"illegal int:"<<buf<<G4endl;
  }
  return 0;
}

// --------------------------------------------------------------------
G4int G4UIparameter::ExpectExponent(const char* str)  // used only by IsDouble()
{
  G4int maxExplength;
  if(IsInt(str, maxExplength = 7) != 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

// --------------------------------------------------------------------
G4int G4UIparameter::IsDouble(
  const char* buf)  // see state diagram for this spec.
{
  const char* p = buf;
  switch(*p)
  {
    case '+':
    case '-':
      ++p;
      if(isdigit(*p) != 0)
      {
        while(isdigit((G4int) (*p)) != 0)
        {
          ++p;
        }
        switch(*p)
        {
          case '\0':
            return 1;  // break;
          case 'E':
          case 'e':
            return ExpectExponent(++p);  // break;
          case '.':
            ++p;
            if(*p == '\0')
            {
              return 1;
            }
            if(*p == 'e' || *p == 'E')
            {
              return ExpectExponent(++p);
            }
            if(isdigit(*p) != 0)
            {
              while(isdigit((G4int) (*p)) != 0)
              {
                ++p;
              }
              if(*p == '\0')
              {
                return 1;
              }
              if(*p == 'e' || *p == 'E')
              {
                return ExpectExponent(++p);
              }
            }
            else
            {
              return 0;
            }
            break;
          default:
            return 0;
        }
      }
      if(*p == '.')
      {
        ++p;
        if(isdigit(*p) != 0)
        {
          while(isdigit((G4int) (*p)) != 0)
          {
            ++p;
          }
          if(*p == '\0')
          {
            return 1;
          }
          if(*p == 'e' || *p == 'E')
          {
            return ExpectExponent(++p);
          }
        }
      }
      break;
    case '.':
      ++p;
      if(isdigit(*p) != 0)
      {
        while(isdigit((G4int) (*p)) != 0)
        {
          ++p;
        }
        if(*p == '\0')
        {
          return 1;
        }
        if(*p == 'e' || *p == 'E')
        {
          return ExpectExponent(++p);
        }
      }
      break;
    default:  // digit is expected
      if(isdigit(*p) != 0)
      {
        while(isdigit((G4int) (*p)) != 0)
        {
          ++p;
        }
        if(*p == '\0')
        {
          return 1;
        }
        if(*p == 'e' || *p == 'E')
        {
          return ExpectExponent(++p);
        }
        if(*p == '.')
        {
          ++p;
          if(*p == '\0')
          {
            return 1;
          }
          if(*p == 'e' || *p == 'E')
          {
            return ExpectExponent(++p);
          }
          if(isdigit(*p) != 0)
          {
            while(isdigit((G4int) (*p)) != 0)
            {
              ++p;
            }
            if(*p == '\0')
            {
              return 1;
            }
            if(*p == 'e' || *p == 'E')
            {
              return ExpectExponent(++p);
            }
          }
        }
      }
  }
  return 0;
}

// ------------------ syntax node functions  ------------------

yystype G4UIparameter::Expression()
{
  yystype result;
#ifdef DEBUG
  G4cerr << " Expression()" << G4endl;
#endif
  result = LogicalORExpression();
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::LogicalORExpression()
{
  yystype result;
  yystype p;
  p = LogicalANDExpression();
  if(token != LOGICALOR)
  {
    return p;
  }
  if(p.type == CONSTSTRING || p.type == IDENTIFIER)
  {
    G4cerr << "Parameter range: illegal type at '||'" << G4endl;
    paramERR = 1;
  }
  result.I = p.I;
  while(token == LOGICALOR)
  {
    token = Yylex();
    p     = LogicalANDExpression();
    if(p.type == CONSTSTRING || p.type == IDENTIFIER)
    {
      G4cerr << "Parameter range: illegal type at '||'" << G4endl;
      paramERR = 1;
    }
    switch(p.type)
    {
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
  if(token != LOGICALAND)
  {
    return p;
  }
  if(p.type == CONSTSTRING || p.type == IDENTIFIER)
  {
    G4cerr << "Parameter range: illegal type at '&&'" << G4endl;
    paramERR = 1;
  }
  result.I = p.I;
  while(token == LOGICALAND)
  {
    token = Yylex();
    p     = EqualityExpression();
    if(p.type == CONSTSTRING || p.type == IDENTIFIER)
    {
      G4cerr << "Parameter range: illegal type at '&&'" << G4endl;
      paramERR = 1;
    }
    switch(p.type)
    {
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
#ifdef DEBUG
  G4cerr << " EqualityExpression()" << G4endl;
#endif
  result = RelationalExpression();
  if(token == EQ || token == NE)
  {
    operat      = token;
    token       = Yylex();
    arg1        = result;
    arg2        = RelationalExpression();
    result.I    = Eval2(arg1, operat, arg2);  // semantic action
    result.type = CONSTINT;
#ifdef DEBUG
    G4cerr << " return code of Eval2(): " << result.I << G4endl;
#endif
  }
  else
  {
    if(result.type != CONSTINT && result.type != CONSTDOUBLE)
    {
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
#ifdef DEBUG
  G4cerr << " RelationalExpression()" << G4endl;
#endif

  arg1 = AdditiveExpression();
  if(token == GT || token == GE || token == LT || token == LE)
  {
    operat      = token;
    token       = Yylex();
    arg2        = AdditiveExpression();
    result.I    = Eval2(arg1, operat, arg2);  // semantic action
    result.type = CONSTINT;
#ifdef DEBUG
    G4cerr << " return  Eval2(): " << G4endl;
#endif
  }
  else
  {
    result = arg1;
  }
#ifdef DEBUG
  G4cerr << " return RelationalExpression()" << G4endl;
#endif
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::AdditiveExpression()
{
  yystype result;
  result = MultiplicativeExpression();
  if(token != '+' && token != '-')
  {
    return result;
  }
  G4cerr << "Parameter range: operator " << (char) token << " is not supported."
         << G4endl;
  paramERR = 1;
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::MultiplicativeExpression()
{
  yystype result;
  result = UnaryExpression();
  if(token != '*' && token != '/' && token != '%')
  {
    return result;
  }
  G4cerr << "Parameter range: operator " << (char) token << " is not supported."
         << G4endl;
  paramERR = 1;
  return result;
}

// --------------------------------------------------------------------
yystype G4UIparameter::UnaryExpression()
{
  yystype result;
  yystype p;
#ifdef DEBUG
  G4cerr << " UnaryExpression" << G4endl;
#endif
  switch(token)
  {
    case '-':
      token = Yylex();
      p     = UnaryExpression();
      if(p.type == CONSTINT)
      {
        result.I    = -p.I;
        result.type = CONSTINT;
      }
      if(p.type == CONSTLONG)
      {
        result.L    = -p.L;
        result.type = CONSTLONG;
      }
      if(p.type == CONSTDOUBLE)
      {
        result.D    = -p.D;
        result.type = CONSTDOUBLE;
      }
      break;
    case '+':
      token  = Yylex();
      result = UnaryExpression();
      break;
    case '!':
      token = Yylex();
      G4cerr << "Parameter range error: "
             << "operator '!' is not supported (sorry)." << G4endl;
      paramERR = 1;
      result   = UnaryExpression();
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
#ifdef DEBUG
  G4cerr << " PrimaryExpression" << G4endl;
#endif
  switch(token)
  {
    case IDENTIFIER:
      result.S    = yylval.S;
      result.type = token;
      token       = Yylex();
      break;
    case CONSTINT:
      result.I    = yylval.I;
      result.type = token;
      token       = Yylex();
      break;
    case CONSTLONG:
      result.L    = yylval.L;
      result.type = token;
      token       = Yylex();
      break;
    case CONSTDOUBLE:
      result.D    = yylval.D;
      result.type = token;
      token       = Yylex();
      break;
    case '(':
      token  = Yylex();
      result = Expression();
      if(token != ')')
      {
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
  if((arg1.type != IDENTIFIER) && (arg2.type != IDENTIFIER))
  {
    G4cerr << parameterName << ": meaningless comparison " << G4int(arg1.type)
           << " " << G4int(arg2.type) << G4endl;
    paramERR = 1;
  }
  char type = (char)std::toupper(parameterType);
  if(arg1.type == IDENTIFIER)
  {
    switch(type)
    {
      case 'I':
        if(arg2.type == CONSTINT)
        {
          return CompareInt(newVal.I, op, arg2.I);
        }
        else
        {
          G4cerr << "integer operand expected for " << parameterRange << '.'
                 << G4endl;
        }
        break;
      case 'L':
        if(arg2.type == CONSTLONG)
        {
          return CompareLong(newVal.L, op, arg2.L);
        }
        else
        {
          G4cerr << "long int operand expected for " << parameterRange << '.'
                 << G4endl;
        }
        break;
      case 'D':
        if(arg2.type == CONSTDOUBLE)
        {
          return CompareDouble(newVal.D, op, arg2.D);
        }
        else if(arg2.type == CONSTINT)
        {  // integral promotion
          return CompareDouble(newVal.D, op, arg2.I);
        }
        else if(arg2.type == CONSTLONG)
        {
          return CompareDouble(newVal.D, op, arg2.L);
        }
        break;
      default:;
    }
  }
  if(arg2.type == IDENTIFIER)
  {
    switch(type)
    {
      case 'I':
        if(arg1.type == CONSTINT)
        {
          return CompareInt(arg1.I, op, newVal.I);
        }
        else
        {
          G4cerr << "integer operand expected for " << parameterRange << '.'
                 << G4endl;
        }
        break;
      case 'L':
        if(arg1.type == CONSTLONG)
        {
          return CompareLong(arg1.L, op, newVal.L);
        }
        else
        {
          G4cerr << "long int operand expected for " << parameterRange << '.'
                 << G4endl;
        }
        break;
      case 'D':
        if(arg1.type == CONSTDOUBLE)
        {
          return CompareDouble(arg1.D, op, newVal.D);
        }
        else if(arg1.type == CONSTINT)
        {  // integral promotion
          return CompareDouble(arg1.I, op, newVal.D);
        }
        else if(arg1.type == CONSTLONG)
        {  // integral promotion
          return CompareDouble(arg1.L, op, newVal.D);
        }
        break;
      default:;
    }
  }
  G4cerr << "no param name is specified at the param range." << G4endl;
  return 0;
}

// --------------------------------------------------------------------
G4int G4UIparameter::CompareInt(G4int arg1, G4int op, G4int arg2)
{
  G4int result = -1;
  G4String opr;
  switch(op)
  {
    case GT:
      result = static_cast<G4int>(arg1 > arg2);
      opr    = ">";
      break;
    case GE:
      result = static_cast<G4int>(arg1 >= arg2);
      opr    = ">=";
      break;
    case LT:
      result = static_cast<G4int>(arg1 < arg2);
      opr    = "<";
      break;
    case LE:
      result = static_cast<G4int>(arg1 <= arg2);
      opr    = "<=";
      break;
    case EQ:
      result = static_cast<G4int>(arg1 == arg2);
      opr    = "==";
      break;
    case NE:
      result = static_cast<G4int>(arg1 != arg2);
      opr    = "!=";
      break;
    default:
      G4cerr << "Parameter range: error at CompareInt" << G4endl;
      paramERR = 1;
  }
#ifdef DEBUG
  G4cerr << "CompareInt " << arg1 << " " << opr << arg2 << " result: " << result
         << G4endl;
#endif
  return result;
}

// --------------------------------------------------------------------
G4int G4UIparameter::CompareLong(G4long arg1, G4int op, G4long arg2)
{
  G4int result = -1;
  G4String opr;
  switch(op)
  {
    case GT:
      result = static_cast<G4int>(arg1 > arg2);
      opr    = ">";
      break;
    case GE:
      result = static_cast<G4int>(arg1 >= arg2);
      opr    = ">=";
      break;
    case LT:
      result = static_cast<G4int>(arg1 < arg2);
      opr    = "<";
      break;
    case LE:
      result = static_cast<G4int>(arg1 <= arg2);
      opr    = "<=";
      break;
    case EQ:
      result = static_cast<G4int>(arg1 == arg2);
      opr    = "==";
      break;
    case NE:
      result = static_cast<G4int>(arg1 != arg2);
      opr    = "!=";
      break;
    default:
      G4cerr << "Parameter range: error at CompareInt" << G4endl;
      paramERR = 1;
  }
#ifdef DEBUG
  G4cerr << "CompareInt " << arg1 << " " << opr << arg2 << " result: " << result
         << G4endl;
#endif
  return result;
}

// --------------------------------------------------------------------
G4int G4UIparameter::CompareDouble(G4double arg1, G4int op, G4double arg2)
{
  G4int result = -1;
  G4String opr;
  switch(op)
  {
    case GT:
      result = static_cast<G4int>(arg1 > arg2);
      opr    = ">";
      break;
    case GE:
      result = static_cast<G4int>(arg1 >= arg2);
      opr    = ">=";
      break;
    case LT:
      result = static_cast<G4int>(arg1 < arg2);
      opr    = "<";
      break;
    case LE:
      result = static_cast<G4int>(arg1 <= arg2);
      opr    = "<=";
      break;
    case EQ:
      result = static_cast<G4int>(arg1 == arg2);
      opr    = "==";
      break;
    case NE:
      result = static_cast<G4int>(arg1 != arg2);
      opr    = "!=";
      break;
    default:
      G4cerr << "Parameter range: error at CompareDouble" << G4endl;
      paramERR = 1;
  }
#ifdef DEBUG
  G4cerr << "CompareDouble " << arg1 << " " << opr << " " << arg2
         << " result: " << result << G4endl;
#endif
  return result;
}

// --------------------- utility functions --------------------------

tokenNum G4UIparameter::Yylex()  // reads input and returns token number KR486
{                                // (returns EOF)
  G4int c;
  G4String buf;

  while((c = G4UIpGetc()) == ' ' || c == '\t' || c == '\n')
  {
    ;
  }
  if(c == EOF)
  {
    return (tokenNum) EOF;  // KR488
  }
  buf = "";
  if((isdigit(c) != 0) || c == '.')
  {  // I or D
    do
    {
      buf += (unsigned char) c;
      c = G4UIpGetc();
    } while(c == '.' || (isdigit(c) != 0) || c == 'e' || c == 'E' || c == '+' ||
            c == '-');
    G4UIpUngetc(c);
    const char* t = buf;
    std::istringstream is(t);
    if(IsInt(buf.data(), 20) != 0)
    {
      is >> yylval.I;
      return CONSTINT;
    }
    else if(IsDouble(buf.data()) != 0)
    {
      is >> yylval.D;
      return CONSTDOUBLE;
    }
    else
    {
      G4cerr << buf << ": numeric format error." << G4endl;
    }
  }
  buf = "";
  if((isalpha(c) != 0) || c == '_')
  {  // IDENTIFIER
    do
    {
      buf += (unsigned char) c;
    } while((c = G4UIpGetc()) != EOF && ((isalnum(c) != 0) || c == '_'));
    G4UIpUngetc(c);
    if(buf == parameterName)
    {
      yylval.S = buf;
      return IDENTIFIER;
    }
    else
    {
      G4cerr << buf << " is not a parameter name." << G4endl;
      paramERR = 1;
    }
  }
  switch(c)
  {
    case '>':
      return (tokenNum) Follow('=', GE, GT);
    case '<':
      return (tokenNum) Follow('=', LE, LT);
    case '=':
      return (tokenNum) Follow('=', EQ, '=');
    case '!':
      return (tokenNum) Follow('=', NE, '!');
    case '|':
      return (tokenNum) Follow('|', LOGICALOR, '|');
    case '&':
      return (tokenNum) Follow('&', LOGICALAND, '&');
    default:
      return (tokenNum) c;
  }
}

// --------------------------------------------------------------------
G4int G4UIparameter::Follow(G4int expect, G4int ifyes, G4int ifno)
{
  G4int c = G4UIpGetc();
  if(c == expect)
  {
    return ifyes;
  }
  G4UIpUngetc(c);
  return ifno;
}

//------------------ low level routines -----------------------------

G4int G4UIparameter::G4UIpGetc()
{  // emulation of getc()
  G4int length = (G4int)parameterRange.length();
  if(bp < length)
  {
    return parameterRange[bp++];
  }
  else
  {
    return EOF;
  }
}

// --------------------------------------------------------------------
G4int G4UIparameter::G4UIpUngetc(G4int c)
{  // emulation of ungetc()
  if(c < 0)
  {
    return -1;
  }
  if(bp > 0 && c == parameterRange[bp - 1])
  {
    --bp;
  }
  else
  {
    G4cerr << "G4UIpUngetc() failed." << G4endl;
    G4cerr << "bp=" << bp << " c=" << c
           << " pR(bp-1)=" << parameterRange[bp - 1] << G4endl;
    paramERR = 1;
    return -1;
  }
  return 0;
}
// *****  end of CheckNewValue() related code  ******
