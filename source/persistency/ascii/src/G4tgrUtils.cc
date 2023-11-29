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
// G4tgrUtils implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include <iomanip>
#include <set>

#include "G4tgrUtils.hh"

#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4tgrParameterMgr.hh"
#include "G4tgrMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
#include "G4UIcommand.hh"

G4ThreadLocal G4tgrEvaluator* G4tgrUtils::theEvaluator = nullptr;

// --------------------------------------------------------------------
G4tgrUtils::G4tgrUtils()
{
  if(theEvaluator == nullptr)
  {
    theEvaluator = new G4tgrEvaluator;
  }
}

// --------------------------------------------------------------------
G4tgrUtils::~G4tgrUtils()
{
  delete theEvaluator;
  theEvaluator = nullptr;
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::IsSeparator(const char ch)
{
  char nonCharacters[7] = { "()+-*/" };
  for(std::size_t ii = 0; ii < 6; ++ii)
  {
    if(ch == nonCharacters[ii])
    {
      return true;
    }
  }
  return false;
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::IsNumber(const G4String& str)
{
  G4int isnum = 1;
  G4int numE  = 0;
  for(G4int ii = 0; ii < (G4int)str.length(); ++ii)
  {
    if(!isdigit(str[ii]) && (str[ii] != '.') && (str[ii] != '-') &&
       (str[ii] != '+'))
    {
      //--- check for E(xponential)
      if(str[ii] == 'E' || str[ii] == 'e')
      {
        if(ii == 0)
        {
          return 0;
        }
        if(numE != 0 || ii == G4int(str.length() - 1))
        {
          isnum = 0;
          break;
        }
        numE++;
      }
      else
      {
        isnum = 0;
        break;
      }
    }
  }
  return isnum;
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::IsInteger(const G4double val, const G4double precision)
{
  if(G4int(val) / val - 1 > precision)
  {
    return false;
  }
  else
  {
    return true;
  }
}

// --------------------------------------------------------------------
void G4tgrUtils::Dump3v(const G4ThreeVector& vec, const char* msg)
{
  G4cout << msg << std::setprecision(8) << vec << std::setprecision(6)
         << G4endl;
}

// --------------------------------------------------------------------
void G4tgrUtils::Dumprm(const G4RotationMatrix& rm, const char* msg)
{
  G4cout << msg << G4endl << " xx=" << rm.xx() << " yx=" << rm.yx()
         << " zx=" << rm.zx() << G4endl << " xy=" << rm.xy()
         << " yy=" << rm.yy() << " zy=" << rm.zy() << G4endl
         << " xz=" << rm.xz() << " yz=" << rm.yz() << " zz=" << rm.zz()
         << G4endl;
}

// --------------------------------------------------------------------
void G4tgrUtils::DumpVS(const std::vector<G4String>& wl, const char* msg,
                        std::ostream& outs)
{
  outs << msg << G4endl;
  for(auto ite = wl.cbegin(); ite != wl.cend(); ++ite)
  {
    outs << *ite << " ";
  }
  outs << G4endl;
}

// --------------------------------------------------------------------
void G4tgrUtils::DumpVS(const std::vector<G4String>& wl, const char* msg)
{
  DumpVS(wl, msg, G4cout);
}

// --------------------------------------------------------------------
G4String G4tgrUtils::SubColon(const G4String& str)
{
  if(str.find(':') != 0)
  {
    G4String ErrMessage = "Trying to subtract leading colon from a word\n" +
                          G4String("that has no leading colon: ") + str;
    G4Exception("G4tgrUtils::SubColon()", "ParseError", FatalException,
                ErrMessage);
  }
  G4String strt = str.substr(1, str.size() - 1);
  return strt;
}

// --------------------------------------------------------------------
G4String G4tgrUtils::GetString(const G4String& str)
{
  //----------- first check if it is parameter
  const char* cstr = str.c_str();
  if(cstr[0] == '$')
  {
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 3)
    {
      G4cout << " G4tgrUtils::GetString() - Substitute parameter: "
             << G4tgrParameterMgr::GetInstance()->FindParameter(
                  str.substr(1, str.size()))
             << G4endl;
    }
#endif
    return G4tgrParameterMgr::GetInstance()->FindParameter(
      str.substr(1, str.size()));
  }
  else
  {
    return str;
  }
}

// --------------------------------------------------------------------
G4double G4tgrUtils::GetDouble(const G4String& str, G4double unitval)
{
  if(!theEvaluator)
  {
    theEvaluator = new G4tgrEvaluator;
  }
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 3)
  {
    G4cout << "G4tgrUtils::GetDouble() - Processing: " << str
           << " default unit " << unitval << G4endl;
  }
#endif
  if(str == "DBL_MAX")
  {
    return DBL_MAX;
  }
  else if(str == "DBL_MIN")
  {
    return DBL_MIN;
  }
  else if(str == "FLT_MAX")
  {
    return FLT_MAX;
  }
  else if(str == "FLT_MIN")
  {
    return FLT_MIN;
  }
  else if(str == "INT_MAX")
  {
    return INT_MAX;
  }
  else if(str == "INT_MIN")
  {
    return INT_MIN;
  }
  //----- Look for arithmetic symbols, (, )
  const char* cstr = str.c_str();
  std::set<G4int> separators;
  separators.insert(-1);
  G4int strlen = G4int(str.length());
  for(G4int ii = 0; ii < strlen; ++ii)
  {
    char cs = cstr[ii];
    if(cs == '*' || cs == '/' || cs == '(' || cs == ')')
    {
      separators.insert(ii);
    }
    else if(cs == '+' || cs == '-')
    {
      // Check if it is not an exponential
      //
      if((ii < 2) || ((cstr[ii - 1] != 'E') && (cstr[ii - 1] != 'e')) ||
	 !IsNumber(G4String(1,cstr[ii - 2])))
      {
        separators.insert(ii);
      }
    }
  }
  separators.insert(strlen);
  std::string strnew;  // build a new word with Parameters
                       // and units substituted by values
  //----- Process words, defined as characters between two separators
  G4int nUnits = 0;
  std::set<G4int>::const_iterator site, site2;
  site  = separators.cbegin();
  site2 = site;
  ++site2;
  for(; site2 != separators.cend(); ++site, ++site2)
  {
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 3)
    {
      G4cout << "   Loop to find word between " << *site << " " << *site2
             << G4endl;
    }
#endif

    if(*site != -1)
    {
      strnew += str.substr(*site, 1);
    }

    G4int wlen = (*site2) - (*site) - 1;  // do not count contiguous separators
    std::string word;
    if(wlen != 0)
    {
      word = str.substr((*site) + 1, (*site2) - (*site) - 1);
    }
    else
    {
      //--- Check combination of separators
      //--- Check number of parentheses
      continue;
    }

#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 3)
    {
      G4cout << "   Processing word: " << word << G4endl;
    }
#endif
    //----------- first check if it is parameter
    const char* cword = word.c_str();
    if(cword[0] == '$')
    {
      G4String parstr = G4tgrParameterMgr::GetInstance()->FindParameter(
        word.substr(1, word.size()));
      if(parstr.substr(0, 1) == "-")
      {
        strnew += "(";
      }
      strnew += parstr;
      if(parstr.substr(0, 1) == "-")
      {
        strnew += ")";
      }
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 3)
      {
        G4cout << " G4tgrutils::GetDouble() - Param found: " << word
               << " in string " << str << " , substituted by " << parstr
               << G4endl;
      }
#endif
    }
    else
    {
      //----- Get if it is a number
      if(IsNumber(word))
      {
        //--- left separator cannot be ')'
        if((*site != -1) && (cstr[*site] == ')'))
        {
          G4String ErrMessage =
            "There cannot be a ')' before a number: " + word +
            " in string: " + str;
          G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                      ErrMessage);
        }
        //--- right separator cannot be '('
        if((*site2 != strlen) && (cstr[*site2] == '('))
        {
          G4String ErrMessage =
            "There cannot be a '(' after a number: " + word +
            " in string: " + str;
          G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                      ErrMessage);
        }
        strnew += word;

        //------ If it is an string, check if it is a unit
      }
      else
      {
        //--- First character cannot be a digit
        if(isdigit(word[0]))
        {
          G4String ErrMessage =
            "String words cannot start with a digit: " + word +
            " in string: " + str;
          G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                      ErrMessage);
        }

        //----- Check if it is a function
        G4bool bWordOK = false;
        if(G4tgrUtils::IsFunction(word))
        {
          //--- It must be followed by '('
          if((*site2 == strlen) || (cstr[*site2] != '('))
          {
            G4String ErrMessage =
              "There must be a '(' after a function: " + word +
              " in string: " + str;
            G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                        ErrMessage);
          }
          strnew += word;
          bWordOK = true;
          //----- Check if it is a unit
        }
        else if(G4tgrUtils::WordIsUnit(word))
        {
          //--- It must be preceded by a *
          if((*site == -1) || ((cstr[*site] != '*') && (cstr[*site] != '/')))
          {
            G4String ErrMess =
              "There must be a '*' before a unit definition: " + word +
              " in string " + str;
            G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                        ErrMess);
          }
          //--- check that it is indeed a CLHEP unit
          if(G4UnitDefinition::GetValueOf(word) != 0.)
          {
            bWordOK = true;
            nUnits++;
            if(nUnits > 1)
            {
              // G4String ErrMess = "There cannot be two unit definitions: "
              //                  + word + " in string " + str;
              // G4Exception("G4tgrUtils::GetDouble()", "ParseError",
              //             FatalException, ErrMess );
            }
            strnew +=
              G4UIcommand::ConvertToString(G4UnitDefinition::GetValueOf(word));
          }
        }
        if(!bWordOK)
        {
          G4String ErrMess = "String word is not a parameter, nor a unit\n" +
                             G4String("definition nor a function: ") + word +
                             G4String(" in string: ") + str;
          G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                      ErrMess);
        }
      }
    }
  }

  G4double val = theEvaluator->evaluate(strnew.c_str());
  if(theEvaluator->status() != HepTool::Evaluator::OK)
  {
    theEvaluator->print_error(theEvaluator->status());
    G4String ErrMessage = "Evaluator error: " + strnew;
    G4Exception("G4tgrUtils::GetDouble()", "ParseError", FatalException,
                ErrMessage);
  }

  if(nUnits == 0)
  {
    val *= unitval;
  }

#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 3)
  {
    G4cout << " G4tgrUtils::GetDouble() - RESULT= " << val << G4endl
           << "   from string: " << str << " converted to: " << strnew.c_str()
           << " with unit val: " << unitval << G4endl;
  }
#endif

  return val;
}

// --------------------------------------------------------------------
G4int G4tgrUtils::GetInt(const G4String& str)
{
  //----- Convert it to a number (it can be a parameter)
  G4double val = GetDouble(str);

  //----- Check it is an integer
  if(!IsInteger(val))
  {
    G4String ErrMessage = G4String("Trying to get the integer from a number") +
                          G4String(" which is not an integer ") + str;
    G4Exception("G4tgrUtils::GetInt()", "ParseError", FatalException,
                ErrMessage);
  }
  return G4int(val);
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::GetBool(const G4String& str)
{
  G4bool val = false;

  //----------- first check that it is a not number
  if((str == "ON") || (str == "TRUE"))
  {
    val = true;
  }
  else if((str == "OFF") || (str == "FALSE"))
  {
    val = false;
  }
  else
  {
    G4String ErrMessage = G4String("Trying to get a float from a string") +
                          G4String(" which is not 'ON'/'OFF'/'TRUE'/'FALSE' ") +
                          str;
    G4Exception("G4tgrUtils::GetBool()", "ParseError", FatalException,
                ErrMessage);
  }

  return val;
}

// --------------------------------------------------------------------
void G4tgrUtils::CheckWLsize(const std::vector<G4String>& wl,
                             unsigned int nWcheck, WLSIZEtype st,
                             const G4String& methodName)
{
  G4String outStr = methodName + G4String(".  Line read with number of words ");
  unsigned int wlsize = (unsigned int)wl.size();

  G4bool isOK = CheckListSize(wlsize, nWcheck, st, outStr);

  if(!isOK)
  {
    G4String chartmp = G4UIcommand::ConvertToString(G4int(nWcheck));
    outStr += chartmp + G4String(" words");
    DumpVS(wl, outStr.c_str());
    G4String ErrMessage =
      " NUMBER OF WORDS: " + G4UIcommand::ConvertToString(G4int(wlsize));
    G4Exception("G4tgrUtils::CheckWLsize()", "ParseError", FatalException,
                ErrMessage);
  }
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::CheckListSize(unsigned int nWreal, unsigned int nWcheck,
                                 WLSIZEtype st, G4String& outStr)
{
  G4bool isOK = true;
  switch(st)
  {
    case WLSIZE_EQ:
      if(nWreal != nWcheck)
      {
        isOK = false;
        outStr += G4String("not equal than ");
      }
      break;
    case WLSIZE_NE:
      if(nWreal == nWcheck)
      {
        isOK = false;
        outStr += G4String("equal than ");
      }
      break;
    case WLSIZE_LE:
      if(nWreal > nWcheck)
      {
        isOK = false;
        outStr += G4String("greater than ");
      }
      break;
    case WLSIZE_LT:
      if(nWreal >= nWcheck)
      {
        isOK = false;
        outStr += G4String("greater or equal than ");
      }
      break;
    case WLSIZE_GE:
      if(nWreal < nWcheck)
      {
        isOK = false;
        outStr += G4String("less than ");
      }
      break;
    case WLSIZE_GT:
      if(nWreal <= nWcheck)
      {
        isOK = false;
        outStr += G4String("less or equal than ");
      }
      break;
    default:
      G4cerr << " ERROR!! - G4tgrUtils::CheckListSize()" << G4endl
             << "           Type of WLSIZE type not found " << st << G4endl;
      break;
  }

  return isOK;
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::WordIsUnit(const G4String& word)
{
  return !IsNumber(word);
  if(word == "mm" || word == "cm" || word == "m" || word == "km" ||
     word == "millimeter" || word == "centimeter" || word == "meter" ||
     word == "kilometer" || word == "parsec" || word == "micrometer" ||
     word == "nanometer" || word == "angstrom" || word == "fermi" ||
     word == "nm" || word == "um" || word == "pc" || word == "radian" ||
     word == "milliradian" || word == "degree" || word == "rad" ||
     word == "mrad" || word == "deg" || word == "ns" || word == "becquerel" ||
     word == "curie")
  {
    return true;
  }
  else
  {
    return false;
  }
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::IsFunction(const G4String& word)
{
  if(word == "sin" || word == "cos" || word == "tan" || word == "asin" ||
     word == "acos" || word == "atan" || word == "atan2" || word == "sinh" ||
     word == "cosh" || word == "tanh" || word == "asinh" || word == "acosh" ||
     word == "atanh" || word == "sqrt" || word == "exp" || word == "log" ||
     word == "log10" || word == "pow")
  {
    return true;
  }
  else
  {
    return false;
  }
}

// --------------------------------------------------------------------
G4RotationMatrix G4tgrUtils::GetRotationFromDirection(G4ThreeVector dir)
{
  G4RotationMatrix rotation;

  if(std::fabs(dir.mag() - 1.) >
     G4GeometryTolerance::GetInstance()->GetSurfaceTolerance())
  {
    G4String WarMessage = "Direction cosines have been normalized to one.\n" +
                          G4String("They were normalized to ") +
                          G4UIcommand::ConvertToString(dir.mag());
    G4Exception("G4tgrUtils::GetRotationFromDirection()", "WrongArgument",
                JustWarning, WarMessage);
    dir /= dir.mag();
  }
  G4double angx = -std::asin(dir.y());

  // There are always two solutions angx, angy and PI-angx,
  // PI+angy, choose first
  //
  G4double angy;
  if(dir.y() == 1.)
  {
    angy = 0.;
  }
  else if(dir.y() == 0.)
  {
    angy = 0.;
  }
  else
  {
    angy = std::asin(dir.x() / std::sqrt(1 - dir.y() * dir.y()));
  }

  // choose between  angy and PI-angy
  if(dir.z() * std::cos(angx) * std::cos(angy) < 0)
  {
    angy = pi - angy;
  }
  rotation.rotateX(angx);
  rotation.rotateY(angy);

  return rotation;
}

// --------------------------------------------------------------------
G4bool G4tgrUtils::AreWordsEquivalent(const G4String& word1,
                                      const G4String& word2)
{
  G4bool bEqual = true;
  std::vector<std::pair<size_t, size_t>> stringPairs;
  // start of substring, number of characters

  //--- Get string limits between asterisks in word1

  std::size_t cStart = 0;
  for(;;)
  {
    size_t cAster = word1.find("*", cStart);
    if(cAster != std::string::npos)
    {
      if(cAster == cStart)
      {
        if(cAster != 0)
        {
          G4Exception("G4tgrUtils::AreWordsEquivalent()",
                      "A word has two asterisks together, please correct it",
                      FatalException, ("Offending word is: " + word1).c_str());
        }
        else
        {
          // word1 == *
          if(word1.size() == 1)
          {
            return true;
          }
        }
      }
      if(cAster != cStart)
      {
        stringPairs.push_back(
          std::pair<size_t, size_t>(cStart, cAster - cStart));
      }
      cStart = cAster + 1;
    }
    else
    {
      if(cStart == 0)
      {
        //--- If there is no asterisk check if they are the same
        return word1 == word2;
      }
      break;
    }
  }

  //---- Add characters after last asterisk as string pair
  if(cStart <= word1.length())
  {
    if(word1.length() != cStart)
    {
      stringPairs.push_back(
        std::pair<size_t, size_t>(cStart, word1.length() - cStart));
    }
  }

  //--- If there are not asterisk, simple comparison
  if(stringPairs.size() == 0)
  {
    if(word1 == word2)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  //--- Find substrings in word2, in same order as in word1
  cStart = 0;
  for(std::size_t ii = 0; ii < stringPairs.size(); ++ii)
  {
    std::pair<size_t, size_t> spair = stringPairs[ii];
    size_t sFound = word2.find(word1.substr(spair.first, spair.second), cStart);
    if(sFound == std::string::npos)
    {
      bEqual = false;
      break;
    }
    else
    {
      //---- If there is no asterisk before first character,
      //     the fisrt string pair found has to start at the first character
      if(spair.first == 0 && sFound != 0)
      {
        bEqual = false;
        break;
        //---- If there is no asterisk after last character,
        //     the last string pair found has to end at the last character
      }
      else if((spair.first + spair.second - 1 == word1.length()) &&
              (sFound + spair.second - 1 != word2.length()))
      {
        bEqual = false;
        break;
      }
      cStart += spair.second;
    }
  }

  return bEqual;
}
