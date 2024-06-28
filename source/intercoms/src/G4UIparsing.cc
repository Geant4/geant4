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
#include "G4UIparsing.hh"

#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UItokenNum.hh"

using namespace G4UItokenNum;

namespace
{
class G4UIRangeChecker
{
  public:
    G4bool DoCheck(const G4UIparameter& p, const char* newValue);
    G4bool DoCheck(const G4UIcommand& cmd, const char* newValue);

  private:
    G4bool RangeCheckImpl(const char* newValue);

    // returns the index of the var name
    unsigned IndexOf(const char* nam);

    // returns 1 or 0
    unsigned IsParameter(const char* nam);

    // --- the following is used by CheckNewValue() -------
    using yystype = G4UItokenNum::yystype;
    using tokenNum = G4UItokenNum::tokenNum;

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
    G4int Follow(G4int expect, G4int ifyes, G4int ifno);

  private:
    // Data from param/cmd
    G4String rangeExpression;
    G4String commandName;
    std::vector<const G4UIparameter*> parameter;

    //------------ CheckNewValue() related data members ---------------
    G4int bp = 0;  // current index in rangeExpression
    tokenNum token = G4UItokenNum::NONE;
    yystype yylval;
    std::vector<yystype> newVal;
    G4int paramERR = 0;
};

// --------------------------------------------------------------------
// INLINE DEFINITIONS
// --------------------------------------------------------------------
inline G4bool G4UIRangeChecker::DoCheck(const G4UIparameter& p, const char* newValue)
{
  // Copy data needed from parameter
  rangeExpression = p.GetParameterRange();

  commandName = p.GetParameterName();
  parameter.resize(1);
  newVal.resize(1);
  parameter[0] = &p;

  return RangeCheckImpl(newValue);
}

inline G4bool G4UIRangeChecker::DoCheck(const G4UIcommand& cmd, const char* newValue)
{
  // Copy data needed from cmd
  rangeExpression = cmd.GetRange();

  commandName = cmd.GetCommandName();
  parameter.resize(cmd.GetParameterEntries());
  newVal.resize(parameter.size());
  for (G4int i = 0; i < (G4int)parameter.size(); ++i) {
    parameter[i] = cmd.GetParameter(i);
  }

  return RangeCheckImpl(newValue);
}

inline G4bool G4UIRangeChecker::RangeCheckImpl(const char* newValue)
{
  if (rangeExpression.empty()) {
    return true;
  }

  yystype result;
  bp = 0;  // reset buffer pointer for G4UIpGetc()
  std::istringstream is(newValue);
  for (unsigned i = 0; i < parameter.size(); ++i) {
    char type = (char)std::toupper(parameter[i]->GetParameterType());
    switch (type) {
      case 'D':
        is >> newVal[i].D;
        break;
      case 'I':
        is >> newVal[i].I;
        break;
      case 'L':
        is >> newVal[i].L;
        break;
      case 'S':
        is >> newVal[i].S;
        break;
      case 'B':
        is >> newVal[i].C;
        break;
      default:;
    }
  }
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

inline unsigned G4UIRangeChecker::IndexOf(const char* nam)
{
  for (unsigned i = 0; i < parameter.size(); ++i) {
    if (parameter[i]->GetParameterName() == nam) {
      return i;
    }
  }
  paramERR = 1;
  G4cerr << "parameter name:" << nam << " not found." << G4endl;
  return 0;
}

inline unsigned G4UIRangeChecker::IsParameter(const char* nam)
{
  for (auto& i : parameter) {
    if (i->GetParameterName() == nam) {
      return 1;
    }
  }
  return 0;
}

// ------------------ syntax node functions  ------------------

inline yystype G4UIRangeChecker::Expression()
{
  return LogicalORExpression();
}

// --------------------------------------------------------------------
inline yystype G4UIRangeChecker::LogicalORExpression()
{
  yystype result;
  yystype p = LogicalANDExpression();
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
inline yystype G4UIRangeChecker::LogicalANDExpression()
{
  yystype result;
  yystype p = EqualityExpression();
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
inline yystype G4UIRangeChecker::EqualityExpression()
{
  yystype arg1, arg2;
  G4int operat;
  yystype result = RelationalExpression();
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
inline yystype G4UIRangeChecker::RelationalExpression()
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
inline yystype G4UIRangeChecker::AdditiveExpression()
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
inline yystype G4UIRangeChecker::MultiplicativeExpression()
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
inline yystype G4UIRangeChecker::UnaryExpression()
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
inline yystype G4UIRangeChecker::PrimaryExpression()
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

//---------------- semantic routines ----------------------------------

inline G4int G4UIRangeChecker::Eval2(const yystype& arg1, G4int op, const yystype& arg2)
{
  if ((arg1.type != IDENTIFIER) && (arg2.type != IDENTIFIER)) {
    G4cerr << commandName << ": meaningless comparison" << G4int(arg1.type) << " "
           << G4int(arg2.type) << G4endl;
    paramERR = 1;
  }

  // Might want to also disallow comparison of same identifiers?

  char newValtype;
  if (arg1.type == IDENTIFIER) {
    unsigned i = IndexOf(arg1.S);
    newValtype = (char)std::toupper(parameter[i]->GetParameterType());
    switch (newValtype) {
      case 'I':
        if (arg2.type == CONSTINT) {
          return G4UIparsing::CompareInt(newVal[i].I, op, arg2.I, paramERR);
          //===================================================================
          // MA - 2018.07.23
        }
        else if (arg2.type == IDENTIFIER) {
          unsigned iii = IndexOf(arg2.S);
          char newValtype2 = (char)std::toupper(parameter[iii]->GetParameterType());
          if (newValtype2 == 'I') {
            return G4UIparsing::CompareInt(newVal[i].I, op, newVal[iii].I, paramERR);
          }
          if (newValtype2 == 'L') {
            G4cerr << "Warning : Integer is compared with long int : " << rangeExpression << G4endl;
            return G4UIparsing::CompareLong(newVal[i].I, op, newVal[iii].L, paramERR);
          }
          if (newValtype2 == 'D') {
            G4cerr << "Warning : Integer is compared with double : " << rangeExpression << G4endl;
            return G4UIparsing::CompareDouble(newVal[i].I, op, newVal[iii].D, paramERR);
          }
          //===================================================================
        }
        else {
          G4cerr << "integer operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'L':
        if (arg2.type == CONSTINT) {
          return G4UIparsing::CompareLong(newVal[i].L, op, arg2.I, paramERR);
        }
        else if (arg2.type == CONSTLONG) {
          return G4UIparsing::CompareLong(newVal[i].L, op, arg2.L, paramERR);
        }
        else if (arg2.type == IDENTIFIER) {
          unsigned iii = IndexOf(arg2.S);
          char newValtype2 = (char)std::toupper(parameter[iii]->GetParameterType());
          if (newValtype2 == 'I') {
            return G4UIparsing::CompareLong(newVal[i].L, op, newVal[iii].I, paramERR);
          }
          if (newValtype2 == 'L') {
            return G4UIparsing::CompareLong(newVal[i].L, op, newVal[iii].L, paramERR);
          }
          if (newValtype2 == 'D') {
            G4cerr << "Warning : Long int is compared with double : " << rangeExpression << G4endl;
            return G4UIparsing::CompareDouble(newVal[i].L, op, newVal[iii].D, paramERR);
          }
          //===================================================================
        }
        else {
          G4cerr << "integer operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'D':
        if (arg2.type == CONSTDOUBLE) {
          return G4UIparsing::CompareDouble(newVal[i].D, op, arg2.D, paramERR);
        }
        else if (arg2.type == CONSTINT) {  // integral promotion
          return G4UIparsing::CompareDouble(newVal[i].D, op, arg2.I, paramERR);
          //===================================================================
          // MA - 2018.07.23
        }
        else if (arg2.type == CONSTLONG) {
          return G4UIparsing::CompareDouble(newVal[i].D, op, arg2.L, paramERR);
        }
        else if (arg2.type == IDENTIFIER) {
          unsigned iii = IndexOf(arg2.S);
          char newValtype2 = (char)std::toupper(parameter[iii]->GetParameterType());
          if (newValtype2 == 'I') {
            return G4UIparsing::CompareDouble(newVal[i].D, op, newVal[iii].I, paramERR);
          }
          if (newValtype2 == 'L') {
            return G4UIparsing::CompareDouble(newVal[i].D, op, newVal[iii].L, paramERR);
          }
          if (newValtype2 == 'D') {
            return G4UIparsing::CompareDouble(newVal[i].D, op, newVal[iii].D, paramERR);
          }
          //===================================================================
        }
        break;
      default:;
    }
  }
  if (arg2.type == IDENTIFIER) {
    unsigned i = IndexOf(arg2.S);
    newValtype = (char)std::toupper(parameter[i]->GetParameterType());
    switch (newValtype) {
      case 'I':
        if (arg1.type == CONSTINT) {
          return G4UIparsing::CompareInt(arg1.I, op, newVal[i].I, paramERR);
        }
        else {
          G4cerr << "integer operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'L':
        if (arg1.type == CONSTLONG) {
          return G4UIparsing::CompareLong(arg1.L, op, newVal[i].L, paramERR);
        }
        else {
          G4cerr << "long int operand expected for " << rangeExpression << '.' << G4endl;
        }
        break;
      case 'D':
        if (arg1.type == CONSTDOUBLE) {
          return G4UIparsing::CompareDouble(arg1.D, op, newVal[i].D, paramERR);
        }
        else if (arg1.type == CONSTINT) {  // integral promotion
          return G4UIparsing::CompareDouble(arg1.I, op, newVal[i].D, paramERR);
        }
        break;
      default:;
    }
  }
  return 0;
}

// --------------------- utility functions ----------------------------

inline tokenNum G4UIRangeChecker::Yylex()  // reads input and returns token number, KR486
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
    if (IsParameter(buf) != 0u) {
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
inline G4int G4UIRangeChecker::Follow(G4int expect, G4int ifyes, G4int ifno)
{
  G4int c = G4UIpGetc();
  if (c == expect) {
    return ifyes;
  }
  G4UIpUngetc(c);
  return ifno;
}

//------------------ low level routines -------------------------------

inline G4int G4UIRangeChecker::G4UIpGetc()
{  // emulation of getc()
  std::size_t length = rangeExpression.length();
  if (bp < (G4int)length) {
    return rangeExpression[bp++];
  }

  return EOF;
}

// --------------------------------------------------------------------
inline G4int G4UIRangeChecker::G4UIpUngetc(G4int c)
{  // emulation of ungetc()
  if (c < 0) {
    return -1;
  }
  if (bp > 0 && c == rangeExpression[bp - 1]) {
    --bp;
  }
  else {
    paramERR = 1;
    return -1;
  }
  return 0;
}
}  // namespace

namespace G4UIparsing
{
G4bool RangeCheck(const G4UIparameter& p, const char* value)
{
  G4UIRangeChecker r;
  return r.DoCheck(p, value);
}

G4bool RangeCheck(const G4UIcommand& p, const char* value)
{
  G4UIRangeChecker r;
  return r.DoCheck(p, value);
}
}  // namespace G4UIRangeCheck
