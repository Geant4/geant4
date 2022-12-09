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
// G4UIcommand
//
// Author: Makoto Asai (SLAC), 1998
// --------------------------------------------------------------------

#include "G4UIcommand.hh"
#include "G4UImessenger.hh"
#include "G4UImanager.hh"
#include "G4UIcommandStatus.hh"
#include "G4StateManager.hh"
#include "G4UnitsTable.hh"
#include "G4Tokenizer.hh"
#include "G4ios.hh"
#include <sstream>
#include <iomanip>

#include "G4UItokenNum.hh"

#include "G4Threading.hh"

using namespace G4UItokenNum;

// --------------------------------------------------------------------
G4UIcommand::G4UIcommand(const char* theCommandPath,
                         G4UImessenger* theMessenger, G4bool tBB)
  : toBeBroadcasted(tBB)
  , messenger(theMessenger)
{
  G4String comStr = theCommandPath;
  G4UIcommandCommonConstructorCode(comStr);
  availabelStateList.clear();
  availabelStateList.push_back(G4State_PreInit);
  availabelStateList.push_back(G4State_Init);
  availabelStateList.push_back(G4State_Idle);
  availabelStateList.push_back(G4State_GeomClosed);
  availabelStateList.push_back(G4State_EventProc);
  availabelStateList.push_back(G4State_Abort);
}

// --------------------------------------------------------------------
void G4UIcommand::G4UIcommandCommonConstructorCode(const char* theCommandPath)
{
  commandPath            = theCommandPath;
  commandName            = theCommandPath;
  G4int commandNameIndex = (G4int)commandName.rfind('/');
  commandName.erase(0, commandNameIndex + 1);
#ifdef G4MULTITHREADED
  if((messenger != nullptr) && messenger->CommandsShouldBeInMaster() &&
     G4Threading::IsWorkerThread())
  {
    toBeBroadcasted = false;
    G4UImanager::GetMasterUIpointer()->AddNewCommand(this);
  }
  else
  {
    G4UImanager::GetUIpointer()->AddNewCommand(this);
  }
#else
  G4UImanager::GetUIpointer()->AddNewCommand(this);
#endif
}

// --------------------------------------------------------------------
void G4UIcommand::SetCommandType(CommandType typ)
{
  if(messenger==nullptr)
  {  // this must be a directory
    if(typ != CmdDirectory)
    {
      G4ExceptionDescription ed;
      ed << "A UI command <" << commandPath
         << "> is defined without vaild messenger.";
      G4Exception("G4UIcommand::SetCommandType","UI2031",
         FatalException,ed);
    }
    else if(commandPath.back() != '/')
    {
      G4ExceptionDescription ed;
      ed << "G4UIcommand Warning : \n" 
         << "  <" << commandPath << "> must be a directory." 
         << "  '/' is appended.";
      G4Exception("G4UIcommand::SetCommandType","UI2032",
         JustWarning,ed);
      commandPath += "/";
    }
  }
  commandType = typ;
}

// --------------------------------------------------------------------
G4UIcommand::~G4UIcommand()
{
  G4UImanager* fUImanager = G4UImanager::GetUIpointer();
  if(fUImanager != nullptr)
  {
    fUImanager->RemoveCommand(this);
  }

  std::size_t n_parameterEntry = parameter.size();
  for(std::size_t i_thParameter = 0; i_thParameter < n_parameterEntry;
      ++i_thParameter)
  {
    delete parameter[i_thParameter];
  }
  parameter.clear();
}

// --------------------------------------------------------------------
G4bool G4UIcommand::operator==(const G4UIcommand& right) const
{
  return (commandPath == right.GetCommandPath());
}

// --------------------------------------------------------------------
G4bool G4UIcommand::operator!=(const G4UIcommand& right) const
{
  return (commandPath != right.GetCommandPath());
}

// --------------------------------------------------------------------
G4int G4UIcommand::DoIt(G4String parameterList)
{
  G4String correctParameters;
  std::size_t n_parameterEntry = parameter.size();
  if(n_parameterEntry != 0)
  {
    G4String aToken;
    G4String correctToken;
    G4Tokenizer parameterToken(parameterList);
    for(std::size_t i_thParameter = 0; i_thParameter < n_parameterEntry;
        ++i_thParameter)
    {
      if(i_thParameter > 0)
      {
        correctParameters.append(" ");
      }
      aToken = parameterToken();
      if(aToken.length() > 0 && aToken[0] == '"')
      {
        while(aToken.back() != '"' ||
              (aToken.length() == 1 && aToken[0] == '"'))
        {
          G4String additionalToken = parameterToken();
          if(additionalToken.empty())
          {
            return G4int(fParameterUnreadable + i_thParameter);
          }
          aToken += " ";
          aToken += additionalToken;
        }
      }
      else if(i_thParameter == n_parameterEntry - 1 &&
              parameter[i_thParameter]->GetParameterType() == 's')
      {
        G4String anotherToken;
        while(!((anotherToken = parameterToken()).empty()))
        {
          std::size_t idxs = anotherToken.find("#");
          if(idxs == std::string::npos)
          {
            aToken += " ";
            aToken += anotherToken;
          }
          else if(idxs > 0)
          {
            aToken += " ";
            aToken += anotherToken.substr(0, idxs);
            break;
          }
          else
          {
            break;
          }
        }
      }

      if(aToken.empty() || aToken == "!")
      {
        if(parameter[i_thParameter]->IsOmittable())
        {
          if(parameter[i_thParameter]->GetCurrentAsDefault())
          {
            G4Tokenizer cvSt(messenger->GetCurrentValue(this));
            G4String parVal;
            for(std::size_t ii = 0; ii < i_thParameter; ++ii)
            {
              parVal = cvSt();
              if(parVal[0] == '"')
              {
                while(parVal.back() != '"')
                {
                  G4String additionalToken = cvSt();
                  if(additionalToken.empty())
                  {
                    return G4int(fParameterUnreadable + i_thParameter);
                  }
                  parVal += " ";
                  parVal += additionalToken;
                }
              }
            }
            G4String aCVToken = cvSt();
            if(aCVToken[0] == '"')
            {
              while(aCVToken.back() != '"')
              {
                G4String additionalToken = cvSt();
                if(additionalToken.empty())
                {
                  return G4int(fParameterUnreadable + i_thParameter);
                }
                aCVToken += " ";
                aCVToken += additionalToken;
              }
            }
            correctParameters.append(aCVToken);
          }
          else
          {
            correctParameters.append(
              parameter[i_thParameter]->GetDefaultValue());
          }
        }
        else
        {
          return G4int(fParameterUnreadable + i_thParameter);
        }
      }
      else
      {
        G4int stat = parameter[i_thParameter]->CheckNewValue(aToken);
        if(stat != 0)
        {
          return stat + G4int(i_thParameter);
        }
        correctParameters.append(aToken);
      }
    }
  }

  if(CheckNewValue(correctParameters) != 0)
  {
    return fParameterOutOfRange + 99;
  }

  if(workerThreadOnly && G4Threading::IsMasterThread())
  {
    return 0;
  }

  messenger->SetNewValue(this, correctParameters);
  return 0;
}

// --------------------------------------------------------------------
G4String G4UIcommand::GetCurrentValue()
{
  return messenger->GetCurrentValue(this);
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1)
{
  availabelStateList.clear();
  availabelStateList.push_back(s1);
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2)
{
  availabelStateList.clear();
  availabelStateList.push_back(s1);
  availabelStateList.push_back(s2);
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2,
                                     G4ApplicationState s3)
{
  availabelStateList.clear();
  availabelStateList.push_back(s1);
  availabelStateList.push_back(s2);
  availabelStateList.push_back(s3);
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2,
                                     G4ApplicationState s3,
                                     G4ApplicationState s4)
{
  availabelStateList.clear();
  availabelStateList.push_back(s1);
  availabelStateList.push_back(s2);
  availabelStateList.push_back(s3);
  availabelStateList.push_back(s4);
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1,
                                     G4ApplicationState s2,
                                     G4ApplicationState s3,
                                     G4ApplicationState s4,
                                     G4ApplicationState s5)
{
  availabelStateList.clear();
  availabelStateList.push_back(s1);
  availabelStateList.push_back(s2);
  availabelStateList.push_back(s3);
  availabelStateList.push_back(s4);
  availabelStateList.push_back(s5);
}

// --------------------------------------------------------------------
G4bool G4UIcommand::IsAvailable()
{
  G4bool av = false;
  G4ApplicationState currentState =
    G4StateManager::GetStateManager()->GetCurrentState();

  std::size_t nState = availabelStateList.size();
  for(std::size_t i = 0; i < nState; ++i)
  {
    if(availabelStateList[i] == currentState)
    {
      av = true;
      break;
    }
  }

  return av;
}

// --------------------------------------------------------------------
G4double G4UIcommand::ValueOf(const char* unitName)
{
  G4double value = 0.;
  value = G4UnitDefinition::GetValueOf(unitName);
  return value;
}

// --------------------------------------------------------------------
G4String G4UIcommand::CategoryOf(const char* unitName)
{
  return G4UnitDefinition::GetCategory(unitName);
}

// --------------------------------------------------------------------
G4String G4UIcommand::UnitsList(const char* unitCategory)
{
  G4String retStr;
  G4UnitsTable& UTbl = G4UnitDefinition::GetUnitsTable();
  std::size_t i;
  for(i = 0; i < UTbl.size(); ++i)
  {
    if(UTbl[i]->GetName() == unitCategory)
    {
      break;
    }
  }
  if(i == UTbl.size())
  {
    G4cerr << "Unit category <" << unitCategory << "> is not defined."
           << G4endl;
    return retStr;
  }
  G4UnitsContainer& UCnt = UTbl[i]->GetUnitsList();
  retStr                 = UCnt[0]->GetSymbol();
  std::size_t je         = UCnt.size();
  for(std::size_t j = 1; j < je; ++j)
  {
    retStr += " ";
    retStr += UCnt[j]->GetSymbol();
  }
  for(std::size_t k = 0; k < je; ++k)
  {
    retStr += " ";
    retStr += UCnt[k]->GetName();
  }
  return retStr;
}

// --------------------------------------------------------------------
void G4UIcommand::List()
{
  G4cout << G4endl;
  G4cout << G4endl;
  if(commandPath.back() != '/')
  {
    G4cout << "Command " << commandPath << G4endl;
  }
  if(workerThreadOnly)
  {
    G4cout << "    ---- available only in worker thread" << G4endl;
  }
  G4cout << "Guidance :" << G4endl;
  std::size_t n_guidanceEntry = commandGuidance.size();
  for(std::size_t i_thGuidance = 0; i_thGuidance < n_guidanceEntry; ++i_thGuidance)
  {
    G4cout << commandGuidance[i_thGuidance] << G4endl;
  }
  if(!rangeString.empty())
  {
    G4cout << " Range of parameters : " << rangeString << G4endl;
  }
  std::size_t n_parameterEntry = parameter.size();
  if(n_parameterEntry > 0)
  {
    for(std::size_t i_thParameter = 0; i_thParameter < n_parameterEntry;
        ++i_thParameter)
    {
      parameter[i_thParameter]->List();
    }
  }
  G4cout << G4endl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4bool boolVal)
{
  G4String vl = "0";
  if(boolVal)
  {
    vl = "1";
  }
  return vl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4int intValue)
{
  std::ostringstream os;
  os << intValue;
  G4String vl = os.str();
  return vl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4long longValue)
{
  std::ostringstream os;
  os << longValue;
  G4String vl = os.str();
  return vl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4double doubleValue)
{
  std::ostringstream os;
  if(G4UImanager::DoublePrecisionStr())
  {
    os << std::setprecision(17) << doubleValue;
  }
  else
  {
    os << doubleValue;
  }
  G4String vl = os.str();
  return vl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4double doubleValue,
                                      const char* unitName)
{
  G4String unt = unitName;
  G4double uv  = ValueOf(unitName);

  std::ostringstream os;
  if(G4UImanager::DoublePrecisionStr())
  {
    os << std::setprecision(17) << doubleValue / uv << " " << unitName;
  }
  else
  {
    os << doubleValue / uv << " " << unitName;
  }
  G4String vl = os.str();
  return vl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(const G4ThreeVector& vec)
{
  std::ostringstream os;
  if(G4UImanager::DoublePrecisionStr())
  {
    os << std::setprecision(17) << vec.x() << " " << vec.y() << " " << vec.z();
  }
  else
  {
    os << vec.x() << " " << vec.y() << " " << vec.z();
  }
  G4String vl = os.str();
  return vl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(const G4ThreeVector& vec,
                                      const char* unitName)
{
  G4String unt = unitName;
  G4double uv  = ValueOf(unitName);

  std::ostringstream os;
  if(G4UImanager::DoublePrecisionStr())
  {
    os << std::setprecision(17) << vec.x() / uv << " " << vec.y() / uv << " "
       << vec.z() / uv << " " << unitName;
  }
  else
  {
    os << vec.x() / uv << " " << vec.y() / uv << " " << vec.z() / uv << " "
       << unitName;
  }
  G4String vl = os.str();
  return vl;
}

// --------------------------------------------------------------------
G4bool G4UIcommand::ConvertToBool(const char* st)
{
  G4String v = G4StrUtil::to_upper_copy(st);
  G4bool vl = false;
  if(v == "Y" || v == "YES" || v == "1" || v == "T" || v == "TRUE")
  {
    vl = true;
  }
  return vl;
}

// --------------------------------------------------------------------
G4int G4UIcommand::ConvertToInt(const char* st)
{
  G4int vl;
  std::istringstream is(st);
  is >> vl;
  return vl;
}

// --------------------------------------------------------------------
G4long G4UIcommand::ConvertToLongInt(const char* st)
{
  G4long vl;
  std::istringstream is(st);
  is >> vl;
  return vl;
}

// --------------------------------------------------------------------
G4double G4UIcommand::ConvertToDouble(const char* st)
{
  G4double vl;
  std::istringstream is(st);
  is >> vl;
  return vl;
}

// --------------------------------------------------------------------
G4double G4UIcommand::ConvertToDimensionedDouble(const char* st)
{
  G4double vl;
  char unts[30];

  std::istringstream is(st);
  is >> vl >> unts;
  G4String unt = unts;

  return (vl * ValueOf(unt));
}

// --------------------------------------------------------------------
G4ThreeVector G4UIcommand::ConvertTo3Vector(const char* st)
{
  G4double vx;
  G4double vy;
  G4double vz;
  std::istringstream is(st);
  is >> vx >> vy >> vz;
  return G4ThreeVector(vx, vy, vz);
}

// --------------------------------------------------------------------
G4ThreeVector G4UIcommand::ConvertToDimensioned3Vector(const char* st)
{
  G4double vx;
  G4double vy;
  G4double vz;
  char unts[30];
  std::istringstream is(st);
  is >> vx >> vy >> vz >> unts;
  G4String unt = unts;
  G4double uv  = ValueOf(unt);
  return G4ThreeVector(vx * uv, vy * uv, vz * uv);
}

// ----- the following is used by CheckNewValue()  --------------------

#include <ctype.h>  // isalpha(), toupper()

G4int G4UIcommand::CheckNewValue(const char* newValue)
{
  yystype result;
  // if( TypeCheck(newValue) == 0 ) return 1;
  if(!rangeString.empty())
  {
    if(RangeCheck(newValue) == 0)
    {
      return fParameterOutOfRange;
    }
  }
  return 0;  // succeeded
}

// --------------------------------------------------------------------
G4int G4UIcommand::TypeCheck(const char* t)
{
  G4String aNewValue;
  char type;
  std::istringstream is(t);
  for(auto& i : parameter)
  {
    is >> aNewValue;
    type = (char)std::toupper(i->GetParameterType());
    switch(type)
    {
      case 'D':
        if(IsDouble(aNewValue) == 0)
        {
          G4cerr << aNewValue << ": double value expected." << G4endl;
          return 0;
        }
        break;
      case 'I':
        if(IsInt(aNewValue, 10) == 0)
        {
          G4cerr << aNewValue << ": integer expected." << G4endl;
          return 0;
        }
        break;
      case 'L':
        if(IsInt(aNewValue, 20) == 0)
        {
          G4cerr << aNewValue << ": long int expected." << G4endl;
          return 0;
        }
        break;
      case 'S':
        break;
      case 'B':
        G4StrUtil::to_upper(aNewValue);
        if(aNewValue == "Y" || aNewValue == "N" || aNewValue == "YES" ||
           aNewValue == "NO" || aNewValue == "1" || aNewValue == "0" ||
           aNewValue == "T" || aNewValue == "F" || aNewValue == "TRUE" ||
           aNewValue == "FALSE")
        {
          return 1;
        }
        else
        {
          return 0;
        }
        break;
      default:;
    }
  }
  return 1;
}

// --------------------------------------------------------------------
G4int G4UIcommand::IsInt(const char* buf, short maxDigits)
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
G4int G4UIcommand::ExpectExponent(const char* str)  // used only by IsDouble()
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
G4int G4UIcommand::IsDouble(const char* buf)
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
            return 1;
            // break;
          case 'E':
          case 'e':
            return ExpectExponent(++p);
            // break;
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

// --------------------------------------------------------------------
G4int G4UIcommand::RangeCheck(const char* t)
{
  yystype result;
  char type;
  bp = 0;  // reset buffer pointer for G4UIpGetc()
  std::istringstream is(t);
  for(unsigned i = 0; i < parameter.size(); ++i)
  {
    type = (char)std::toupper(parameter[i]->GetParameterType());
    switch(type)
    {
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
  G4cerr << "parameter out of range: " << rangeString << G4endl;
  return 0;
}

// ------------------ syntax node functions  ------------------

yystype G4UIcommand::Expression()
{
  yystype result;
#ifdef DEBUG
  G4cerr << " Expression()" << G4endl;
#endif
  result = LogicalORExpression();
  return result;
}

// --------------------------------------------------------------------
yystype G4UIcommand::LogicalORExpression()
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
yystype G4UIcommand::LogicalANDExpression()
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
yystype G4UIcommand::EqualityExpression()
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
yystype G4UIcommand::RelationalExpression()
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
    G4cerr << " return code of Eval2(): " << result.I << G4endl;
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
yystype G4UIcommand::AdditiveExpression()
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
yystype G4UIcommand::MultiplicativeExpression()
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
yystype G4UIcommand::UnaryExpression()
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
yystype G4UIcommand::PrimaryExpression()
{
  yystype result;
#ifdef DEBUG
  G4cerr << " primary_exp" << G4endl;
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

//---------------- semantic routines ----------------------------------

G4int G4UIcommand::Eval2(const yystype& arg1, G4int op, const yystype& arg2)
{
  char newValtype;
  if((arg1.type != IDENTIFIER) && (arg2.type != IDENTIFIER))
  {
    G4cerr << commandName << ": meaningless comparison" << G4endl;
    paramERR = 1;
  }

  if(arg1.type == IDENTIFIER)
  {
    unsigned i = IndexOf(arg1.S);
    newValtype = (char)std::toupper(parameter[i]->GetParameterType());
    switch(newValtype)
    {
      case 'I':
        if(arg2.type == CONSTINT)
        {
          return CompareInt(newVal[i].I, op, arg2.I);
          //===================================================================
          // MA - 2018.07.23
        }
        else if(arg2.type == IDENTIFIER)
        {
          unsigned iii     = IndexOf(arg2.S);
          char newValtype2 = (char)std::toupper(parameter[iii]->GetParameterType());
          if(newValtype2 == 'I')
          {
            return CompareInt(newVal[i].I, op, newVal[iii].I);
          }
          else if(newValtype2 == 'L')
          {
            G4cerr << "Warning : Integer is compared with long int : "
                   << rangeString << G4endl;
            return CompareLong(newVal[i].I, op, newVal[iii].L);
          }
          else if(newValtype2 == 'D')
          {
            G4cerr << "Warning : Integer is compared with double : "
                   << rangeString << G4endl;
            return CompareDouble(newVal[i].I, op, newVal[iii].D);
          }
          //===================================================================
        }
        else
        {
          G4cerr << "integer operand expected for " << rangeString << '.'
                 << G4endl;
        }
        break;
      case 'L':
        if(arg2.type == CONSTINT)
        {
          return CompareLong(newVal[i].L, op, arg2.I);
        }
        else if(arg2.type == CONSTLONG)
        {
          return CompareLong(newVal[i].L, op, arg2.L);
        }
        else if(arg2.type == IDENTIFIER)
        {
          unsigned iii     = IndexOf(arg2.S);
          char newValtype2 = (char)std::toupper(parameter[iii]->GetParameterType());
          if(newValtype2 == 'I')
          {
            return CompareLong(newVal[i].L, op, newVal[iii].I);
          }
          if(newValtype2 == 'L')
          {
            return CompareLong(newVal[i].L, op, newVal[iii].L);
          }
          else if(newValtype2 == 'D')
          {
            G4cerr << "Warning : Long int is compared with double : "
                   << rangeString << G4endl;
            return CompareDouble(newVal[i].L, op, newVal[iii].D);
          }
          //===================================================================
        }
        else
        {
          G4cerr << "integer operand expected for " << rangeString << '.'
                 << G4endl;
        }
        break;
      case 'D':
        if(arg2.type == CONSTDOUBLE)
        {
          return CompareDouble(newVal[i].D, op, arg2.D);
        }
        else if(arg2.type == CONSTINT)
        {  // integral promotion
          return CompareDouble(newVal[i].D, op, arg2.I);
          //===================================================================
          // MA - 2018.07.23
        }
        else if(arg2.type == CONSTLONG)
        {
          return CompareDouble(newVal[i].D, op, arg2.L);
        }
        else if(arg2.type == IDENTIFIER)
        {
          unsigned iii     = IndexOf(arg2.S);
          char newValtype2 = (char)std::toupper(parameter[iii]->GetParameterType());
          if(newValtype2 == 'I')
          {
            return CompareDouble(newVal[i].D, op, newVal[iii].I);
          }
          else if(newValtype2 == 'L')
          {
            return CompareDouble(newVal[i].D, op, newVal[iii].L);
          }
          else if(newValtype2 == 'D')
          {
            return CompareDouble(newVal[i].D, op, newVal[iii].D);
          }
          //===================================================================
        }
        break;
      default:;
    }
  }
  if(arg2.type == IDENTIFIER)
  {
    unsigned i = IndexOf(arg2.S);
    newValtype = (char)std::toupper(parameter[i]->GetParameterType());
    switch(newValtype)
    {
      case 'I':
        if(arg1.type == CONSTINT)
        {
          return CompareInt(arg1.I, op, newVal[i].I);
        }
        else
        {
          G4cerr << "integer operand expected for " << rangeString << '.'
                 << G4endl;
        }
        break;
      case 'L':
        if(arg1.type == CONSTLONG)
        {
          return CompareLong(arg1.L, op, newVal[i].L);
        }
        else
        {
          G4cerr << "long int operand expected for " << rangeString << '.'
                 << G4endl;
        }
        break;
      case 'D':
        if(arg1.type == CONSTDOUBLE)
        {
          return CompareDouble(arg1.D, op, newVal[i].D);
        }
        else if(arg1.type == CONSTINT)
        {  // integral promotion
          return CompareDouble(arg1.I, op, newVal[i].D);
        }
        break;
      default:;
    }
  }
  return 0;
}

// --------------------------------------------------------------------
G4int G4UIcommand::CompareInt(G4int arg1, G4int op, G4int arg2)
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
G4int G4UIcommand::CompareLong(G4long arg1, G4int op, G4long arg2)
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
G4int G4UIcommand::CompareDouble(G4double arg1, G4int op, G4double arg2)
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

// --------------------------------------------------------------------
unsigned G4UIcommand::IndexOf(const char* nam)
{
  unsigned i;
  G4String pname;
  for(i = 0; i < parameter.size(); ++i)
  {
    pname = parameter[i]->GetParameterName();
    if(pname == nam)
    {
      return i;
    }
  }
  paramERR = 1;
  G4cerr << "parameter name:" << nam << " not found." << G4endl;
  return 0;
}

// --------------------------------------------------------------------
unsigned G4UIcommand::IsParameter(const char* nam)
{
  G4String pname;
  for(auto& i : parameter)
  {
    pname = i->GetParameterName();
    if(pname == nam)
    {
      return 1;
    }
  }
  return 0;
}

// --------------------- utility functions ----------------------------

tokenNum G4UIcommand::Yylex()  // reads input and returns token number, KR486
{                              // (returns EOF)
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
      buf += (unsigned char)c;
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
    if(IsParameter(buf) != 0u)
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
G4int G4UIcommand::Follow(G4int expect, G4int ifyes, G4int ifno)
{
  G4int c = G4UIpGetc();
  if(c == expect)
  {
    return ifyes;
  }
  G4UIpUngetc(c);
  return ifno;
}

//------------------ low level routines -------------------------------

G4int G4UIcommand::G4UIpGetc()
{  // emulation of getc()
  std::size_t length = rangeString.length();
  if(bp < (G4int)length)
  {
    return rangeString[bp++];
  }
  else
  {
    return EOF;
  }
}

// --------------------------------------------------------------------
G4int G4UIcommand::G4UIpUngetc(G4int c)
{  // emulation of ungetc()
  if(c < 0)
  {
    return -1;
  }
  if(bp > 0 && c == rangeString[bp - 1])
  {
    --bp;
  }
  else
  {
    G4cerr << "G4UIpUngetc() failed." << G4endl;
    G4cerr << "bp=" << bp << " c=" << c << " pR(bp-1)=" << rangeString[bp - 1]
           << G4endl;
    paramERR = 1;
    return -1;
  }
  return 0;
}
