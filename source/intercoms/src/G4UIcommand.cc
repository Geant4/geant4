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

#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4Tokenizer.hh"
#include "G4UIcommandStatus.hh"
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4UIparsing.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

// --------------------------------------------------------------------
G4UIcommand::G4UIcommand(const char* theCommandPath, G4UImessenger* theMessenger, G4bool tBB)
  : toBeBroadcasted(tBB), messenger(theMessenger)
{
  G4String comStr = theCommandPath;
  G4UIcommandCommonConstructorCode(comStr);
  availabelStateList = {G4State_PreInit,    G4State_Init,      G4State_Idle,
                        G4State_GeomClosed, G4State_EventProc, G4State_Abort};
}

// --------------------------------------------------------------------
void G4UIcommand::G4UIcommandCommonConstructorCode(const char* theCommandPath)
{
  commandPath = theCommandPath;
  commandName = theCommandPath;
  auto commandNameIndex = (G4int)commandName.rfind('/');
  commandName.erase(0, commandNameIndex + 1);
#ifdef G4MULTITHREADED
  if ((messenger != nullptr) && messenger->CommandsShouldBeInMaster()
      && G4Threading::IsWorkerThread())
  {
    toBeBroadcasted = false;
    G4UImanager::GetMasterUIpointer()->AddNewCommand(this);
  }
  else {
    G4UImanager::GetUIpointer()->AddNewCommand(this);
  }
#else
  G4UImanager::GetUIpointer()->AddNewCommand(this);
#endif
}

// --------------------------------------------------------------------
void G4UIcommand::SetCommandType(CommandType typ)
{
  if (messenger == nullptr) {  // this must be a directory
    if (typ != CmdDirectory) {
      G4ExceptionDescription ed;
      ed << "A UI command <" << commandPath << "> is defined without vaild messenger.";
      G4Exception("G4UIcommand::SetCommandType", "UI2031", FatalException, ed);
    }
    else if (commandPath.back() != '/') {
      G4ExceptionDescription ed;
      ed << "G4UIcommand Warning : \n"
         << "  <" << commandPath << "> must be a directory."
         << "  '/' is appended.";
      G4Exception("G4UIcommand::SetCommandType", "UI2032", JustWarning, ed);
      commandPath += "/";
    }
  }
  commandType = typ;
}

// --------------------------------------------------------------------
G4UIcommand::~G4UIcommand()
{
  G4UImanager* fUImanager = G4UImanager::GetUIpointer();
  if (fUImanager != nullptr) {
    fUImanager->RemoveCommand(this);
  }

  for (const auto& p : parameter) {
    delete p;
  }
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
  if (n_parameterEntry != 0) {
    G4String aToken;
    G4String correctToken;
    G4Tokenizer parameterToken(parameterList);
    for (std::size_t i_thParameter = 0; i_thParameter < n_parameterEntry; ++i_thParameter) {
      if (i_thParameter > 0) {
        correctParameters.append(" ");
      }
      aToken = parameterToken();
      if (aToken.length() > 0 && aToken[0] == '"') {
        while (aToken.back() != '"' || (aToken.length() == 1 && aToken[0] == '"')) {
          G4String additionalToken = parameterToken();
          if (additionalToken.empty()) {
            return G4int(fParameterUnreadable + i_thParameter);
          }
          aToken += " ";
          aToken += additionalToken;
        }
      }
      else if (i_thParameter == n_parameterEntry - 1
               && parameter[i_thParameter]->GetParameterType() == 's')
      {
        G4String anotherToken;
        while (!((anotherToken = parameterToken()).empty())) {
          std::size_t idxs = anotherToken.find('#');
          if (idxs == std::string::npos) {
            aToken += " ";
            aToken += anotherToken;
          }
          else if (idxs > 0) {
            aToken += " ";
            aToken += anotherToken.substr(0, idxs);
            break;
          }
          else {
            break;
          }
        }
      }

      if (aToken.empty() || aToken == "!") {
        if (parameter[i_thParameter]->IsOmittable()) {
          if (parameter[i_thParameter]->GetCurrentAsDefault()) {
            G4Tokenizer cvSt(messenger->GetCurrentValue(this));
            G4String parVal;
            for (std::size_t ii = 0; ii < i_thParameter; ++ii) {
              parVal = cvSt();
              if (parVal[0] == '"') {
                while (parVal.back() != '"') {
                  G4String additionalToken = cvSt();
                  if (additionalToken.empty()) {
                    return G4int(fParameterUnreadable + i_thParameter);
                  }
                  parVal += " ";
                  parVal += additionalToken;
                }
              }
            }
            G4String aCVToken = cvSt();
            if (aCVToken[0] == '"') {
              while (aCVToken.back() != '"') {
                G4String additionalToken = cvSt();
                if (additionalToken.empty()) {
                  return G4int(fParameterUnreadable + i_thParameter);
                }
                aCVToken += " ";
                aCVToken += additionalToken;
              }
            }
            correctParameters.append(aCVToken);
          }
          else {
            correctParameters.append(parameter[i_thParameter]->GetDefaultValue());
          }
        }
        else {
          return G4int(fParameterUnreadable + i_thParameter);
        }
      }
      else {
        G4int stat = parameter[i_thParameter]->CheckNewValue(aToken);
        if (stat != 0) {
          return stat + G4int(i_thParameter);
        }
        correctParameters.append(aToken);
      }
    }
  }

  if (CheckNewValue(correctParameters) != 0) {
    return fParameterOutOfRange + 99;
  }

  if (workerThreadOnly && G4Threading::IsMasterThread()) {
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
  availabelStateList = {s1};
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1, G4ApplicationState s2)
{
  availabelStateList = {s1, s2};
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1, G4ApplicationState s2,
                                     G4ApplicationState s3)
{
  availabelStateList = {s1, s2, s3};
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1, G4ApplicationState s2,
                                     G4ApplicationState s3, G4ApplicationState s4)
{
  availabelStateList = {s1, s2, s3, s4};
}

// --------------------------------------------------------------------
void G4UIcommand::AvailableForStates(G4ApplicationState s1, G4ApplicationState s2,
                                     G4ApplicationState s3, G4ApplicationState s4,
                                     G4ApplicationState s5)
{
  availabelStateList = {s1, s2, s3, s4, s5};
}

// --------------------------------------------------------------------
G4bool G4UIcommand::IsAvailable()
{
  G4ApplicationState currentState = G4StateManager::GetStateManager()->GetCurrentState();

  for (const auto& s : availabelStateList) {
    if (s == currentState) {
      return true;
    }
  }

  return false;
}

// --------------------------------------------------------------------
G4double G4UIcommand::ValueOf(const char* unitName)
{
  return G4UnitDefinition::GetValueOf(unitName);
}

// --------------------------------------------------------------------
G4String G4UIcommand::CategoryOf(const char* unitName)
{
  return G4UnitDefinition::GetCategory(unitName);
}

// --------------------------------------------------------------------
G4String G4UIcommand::UnitsList(const char* unitCategory)
{
  G4UnitsTable& UTbl = G4UnitDefinition::GetUnitsTable();

  auto ucatIter = std::find_if(std::cbegin(UTbl), std::cend(UTbl), [&unitCategory](const auto& ud) {
    return ud->GetName() == unitCategory;
  });

  if (ucatIter == std::cend(UTbl)) {
    G4cerr << "Unit category <" << unitCategory << "> is not defined." << G4endl;
    return G4String();
  }

  G4String symList;
  G4String nameList;
  G4UnitsContainer& UCnt = (*ucatIter)->GetUnitsList();

  for (const auto& uDef : UCnt) {
    symList += uDef->GetSymbol();
    symList += " ";
    nameList += uDef->GetName();
    nameList += " ";
  }
  symList += nameList;
  G4StrUtil::rstrip(symList);
  return symList;
}

// --------------------------------------------------------------------
void G4UIcommand::List()
{
  G4cout << G4endl;
  G4cout << G4endl;
  if (commandPath.back() != '/') {
    G4cout << "Command " << commandPath << G4endl;
  }
  if (workerThreadOnly) {
    G4cout << "    ---- available only in worker thread" << G4endl;
  }

  G4cout << "Guidance :" << G4endl;
  for (const auto& i_thGuidance : commandGuidance) {
    G4cout << i_thGuidance << G4endl;
  }

  if (!rangeExpression.empty()) {
    G4cout << " Range of parameters : " << rangeExpression << G4endl;
  }

  for (const auto& i_thParameter : parameter) {
    i_thParameter->List();
  }
  G4cout << G4endl;
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4bool boolVal)
{
  return boolVal ? "1" : "0";
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4int intValue)
{
  return G4UIparsing::TtoS(intValue);
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4long longValue)
{
  return G4UIparsing::TtoS(longValue);
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4double doubleValue)
{
  std::ostringstream os;
  if (G4UImanager::DoublePrecisionStr()) {
    os << std::setprecision(17);
  }
  os << doubleValue;
  return os.str();
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(G4double doubleValue, const char* unitName)
{
  std::ostringstream os;
  if (G4UImanager::DoublePrecisionStr()) {
    os << std::setprecision(17);
  }
  os << doubleValue / ValueOf(unitName) << " " << unitName;
  return os.str();
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(const G4ThreeVector& vec)
{
  std::ostringstream os;
  if (G4UImanager::DoublePrecisionStr()) {
    os << std::setprecision(17);
  }
  os << vec.x() << " " << vec.y() << " " << vec.z();
  return os.str();
}

// --------------------------------------------------------------------
G4String G4UIcommand::ConvertToString(const G4ThreeVector& vec, const char* unitName)
{
  G4double uv = ValueOf(unitName);

  std::ostringstream os;
  if (G4UImanager::DoublePrecisionStr()) {
    os << std::setprecision(17);
  }
  os << vec.x() / uv << " " << vec.y() / uv << " " << vec.z() / uv << " " << unitName;
  return os.str();
}

// --------------------------------------------------------------------
G4bool G4UIcommand::ConvertToBool(const char* st)
{
  G4String v = G4StrUtil::to_upper_copy(st);
  return (v == "Y" || v == "YES" || v == "1" || v == "T" || v == "TRUE");
}

// --------------------------------------------------------------------
G4int G4UIcommand::ConvertToInt(const char* st)
{
  return G4UIparsing::StoT<G4int>(st);
}

// --------------------------------------------------------------------
G4long G4UIcommand::ConvertToLongInt(const char* st)
{
  return G4UIparsing::StoT<G4long>(st);
}

// --------------------------------------------------------------------
G4double G4UIcommand::ConvertToDouble(const char* st)
{
  return G4UIparsing::StoT<G4double>(st);
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
  G4double uv = ValueOf(unt);
  return G4ThreeVector(vx * uv, vy * uv, vz * uv);
}

G4int G4UIcommand::CheckNewValue(const char* newValue)
{
  if (!G4UIparsing::RangeCheck(*this, newValue)) {
    return fParameterOutOfRange;
  }
  return 0;  // succeeded
}
