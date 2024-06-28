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
  if (!G4UIparsing::RangeCheck(*this, newValue)) {
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
