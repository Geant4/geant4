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
// G4UIparsing
//
// Utilities for parsing/checking inputs in G4UIparameter/G4UIcommand

#ifndef G4UIparsing_hh
#define G4UIparsing_hh 1

#include "G4UItokenNum.hh"
#include "globals.hh"

#include <cctype>

namespace G4UIparsing
{
// Convert G4String to value of type T
template<typename T>
inline T StoT(const G4String& s)
{
  T vl;
  std::istringstream is(s);
  is >> vl;
  return vl;
}

// Convert value of type T to G4String
template<typename T>
inline G4String TtoS(T value)
{
  std::ostringstream os;
  os << value;
  return os.str();
}

// Return true if `str` parses to an integral number no more than `maxDigit` digits
inline G4bool IsInt(const char* str, short maxDigits)
{
  const char* p = str;
  G4int length = 0;
  if (*p == '+' || *p == '-') {
    ++p;
  }
  if (isdigit((G4int)(*p)) != 0) {
    while (isdigit((G4int)(*p)) != 0) {
      ++p;
      ++length;
    }
    if (*p == '\0') {
      if (length > maxDigits) {
        G4cerr << "digit length exceeds" << G4endl;
        return false;
      }
      return true;
    }
  }
  return false;
}

// Return true if `str` parses to an exponent
//
// A valid exponent is an integer of no more than 7 digits
inline G4bool ExpectExponent(const char* str)
{
  return IsInt(str, 7);
}

inline G4bool IsDouble(const char* str)
{
  const char* p = str;
  switch (*p) {
    case '+':
    case '-':
      ++p;
      if (isdigit(*p) != 0) {
        while (isdigit((G4int)(*p)) != 0) {
          ++p;
        }
        switch (*p) {
          case '\0':
            return true;  // break;
          case 'E':
          case 'e':
            return ExpectExponent(++p);  // break;
          case '.':
            ++p;
            if (*p == '\0') {
              return true;
            }
            if (*p == 'e' || *p == 'E') {
              return ExpectExponent(++p);
            }
            if (isdigit(*p) != 0) {
              while (isdigit((G4int)(*p)) != 0) {
                ++p;
              }
              if (*p == '\0') {
                return true;
              }
              if (*p == 'e' || *p == 'E') {
                return ExpectExponent(++p);
              }
            }
            else {
              return false;
            }
            break;
          default:
            return false;
        }
      }
      if (*p == '.') {
        ++p;
        if (isdigit(*p) != 0) {
          while (isdigit((G4int)(*p)) != 0) {
            ++p;
          }
          if (*p == '\0') {
            return true;
          }
          if (*p == 'e' || *p == 'E') {
            return ExpectExponent(++p);
          }
        }
      }
      break;
    case '.':
      ++p;
      if (isdigit(*p) != 0) {
        while (isdigit((G4int)(*p)) != 0) {
          ++p;
        }
        if (*p == '\0') {
          return true;
        }
        if (*p == 'e' || *p == 'E') {
          return ExpectExponent(++p);
        }
      }
      break;
    default:  // digit is expected
      if (isdigit(*p) != 0) {
        while (isdigit((G4int)(*p)) != 0) {
          ++p;
        }
        if (*p == '\0') {
          return true;
        }
        if (*p == 'e' || *p == 'E') {
          return ExpectExponent(++p);
        }
        if (*p == '.') {
          ++p;
          if (*p == '\0') {
            return true;
          }
          if (*p == 'e' || *p == 'E') {
            return ExpectExponent(++p);
          }
          if (isdigit(*p) != 0) {
            while (isdigit((G4int)(*p)) != 0) {
              ++p;
            }
            if (*p == '\0') {
              return true;
            }
            if (*p == 'e' || *p == 'E') {
              return ExpectExponent(++p);
            }
          }
        }
      }
  }
  return false;
}

// --------------------------------------------------------------------
inline G4int CompareInt(G4int arg1, G4int op, G4int arg2, G4int& errCode)
{
  G4int result = -1;
  switch (op) {
    case G4UItokenNum::GT:
      result = static_cast<G4int>(arg1 > arg2);
      break;
    case G4UItokenNum::GE:
      result = static_cast<G4int>(arg1 >= arg2);
      break;
    case G4UItokenNum::LT:
      result = static_cast<G4int>(arg1 < arg2);
      break;
    case G4UItokenNum::LE:
      result = static_cast<G4int>(arg1 <= arg2);
      break;
    case G4UItokenNum::EQ:
      result = static_cast<G4int>(arg1 == arg2);
      break;
    case G4UItokenNum::NE:
      result = static_cast<G4int>(arg1 != arg2);
      break;
    default:
      G4cerr << "Parameter range: error at CompareInt" << G4endl;
      errCode = 1;
  }
  return result;
}

// --------------------------------------------------------------------
inline G4int CompareLong(G4long arg1, G4int op, G4long arg2, G4int& errCode)
{
  G4int result = -1;
  switch (op) {
    case G4UItokenNum::GT:
      result = static_cast<G4int>(arg1 > arg2);
      break;
    case G4UItokenNum::GE:
      result = static_cast<G4int>(arg1 >= arg2);
      break;
    case G4UItokenNum::LT:
      result = static_cast<G4int>(arg1 < arg2);
      break;
    case G4UItokenNum::LE:
      result = static_cast<G4int>(arg1 <= arg2);
      break;
    case G4UItokenNum::EQ:
      result = static_cast<G4int>(arg1 == arg2);
      break;
    case G4UItokenNum::NE:
      result = static_cast<G4int>(arg1 != arg2);
      break;
    default:
      G4cerr << "Parameter range: error at CompareInt" << G4endl;
      errCode = 1;
  }
  return result;
}

// --------------------------------------------------------------------
inline G4int CompareDouble(G4double arg1, G4int op, G4double arg2, G4int& errCode)
{
  G4int result = -1;
  switch (op) {
    case G4UItokenNum::GT:
      result = static_cast<G4int>(arg1 > arg2);
      break;
    case G4UItokenNum::GE:
      result = static_cast<G4int>(arg1 >= arg2);
      break;
    case G4UItokenNum::LT:
      result = static_cast<G4int>(arg1 < arg2);
      break;
    case G4UItokenNum::LE:
      result = static_cast<G4int>(arg1 <= arg2);
      break;
    case G4UItokenNum::EQ:
      result = static_cast<G4int>(arg1 == arg2);
      break;
    case G4UItokenNum::NE:
      result = static_cast<G4int>(arg1 != arg2);
      break;
    default:
      G4cerr << "Parameter range: error at CompareDouble" << G4endl;
      errCode = 1;
  }
  return result;
}

}  // namespace G4UIparsing

#endif
