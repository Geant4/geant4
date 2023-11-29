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
// G4tgrFileIn implementation
//
// Author: P.Arce, CIEMAT (November 2007)
// --------------------------------------------------------------------

#include "globals.hh"

#include <iostream>
#include <fstream>
#include <sstream>

#include "G4tgrFileIn.hh"
#include "G4tgrMessenger.hh"
#include "G4tgrUtils.hh"
#include "G4UIcommand.hh"

G4ThreadLocal std::vector<G4tgrFileIn*>* G4tgrFileIn::theInstances = nullptr;

// --------------------------------------------------------------------
G4tgrFileIn::G4tgrFileIn()
{
  if(theInstances == nullptr)
  {
    theInstances = new std::vector<G4tgrFileIn*>;
  }
}

// --------------------------------------------------------------------
G4tgrFileIn::~G4tgrFileIn()
{
  delete theInstances;
  theInstances = nullptr;
  /*
    for( auto vfcite = theInstances->cbegin();
              vfcite != theInstances->cend(); ++vfcite)
    {
      delete *vfcite;
    }
  */
}

// --------------------------------------------------------------------
G4tgrFileIn& G4tgrFileIn::GetInstance(const G4String& filename)
{
  if(theInstances == nullptr)
  {
    theInstances = new std::vector<G4tgrFileIn*>;
  }

  std::vector<G4tgrFileIn*>::const_iterator vfcite;
  for(vfcite = theInstances->cbegin(); vfcite != theInstances->cend(); ++vfcite)
  {
    if((*vfcite)->GetName() == filename)
    {
      return *(*vfcite);
    }
  }

  G4tgrFileIn* instance = nullptr;
  if(vfcite == theInstances->cend())
  {
    instance = new G4tgrFileIn(filename);

    instance->theCurrentFile = -1;
    instance->OpenNewFile(filename.c_str());

    theInstances->push_back(instance);
  }

  return *instance;
}

// --------------------------------------------------------------------
void G4tgrFileIn::OpenNewFile(const char* filename)
{
  ++theCurrentFile;
  std::ifstream* fin = new std::ifstream(filename);
  theFiles.push_back(fin);

  theLineNo.push_back(0);

  theNames.push_back(filename);

#ifndef OS_SUN_4_2
  if(!fin->is_open())
  {
    G4String ErrMessage = "Input file does not exist: " + G4String(filename);
    G4Exception("G4tgrFileIn::OpenNewFile()", "InvalidInput", FatalException,
                ErrMessage);
  }
#endif
}

// --------------------------------------------------------------------
G4tgrFileIn& G4tgrFileIn::GetInstanceOpened(const G4String& filename)
{
  G4tgrFileIn& filein = G4tgrFileIn::GetInstance(filename);
  if(filein.GetName() != filename)
  {
    G4String ErrMessage = "File not opened yet: " + filename;
    G4Exception("G4tgrFileIn::GetInstanceOpened()", "InvalidInput",
                FatalException, ErrMessage);
  }
  else
  {
    return filein;
  }
  return filein;  // to avoid compilation warnings
}

// --------------------------------------------------------------------
G4int G4tgrFileIn::GetWordsInLine(std::vector<G4String>& wordlist)
{
  G4int isok = 1;

  //---------- Read a line of file:
  //           NOTE: cannot be read with a istream_iterator,
  //           because it uses G4cout, and then doesn't read '\n'
  //----- Clear wordlist
  G4int wsiz = (G4int)wordlist.size();
  G4int ii;
  for(ii = 0; ii < wsiz; ++ii)
  {
    wordlist.pop_back();
  }

  //---------- Loop lines while there is an ending '\' or line is blank
  const G4int NMAXLIN = 1000;
  char ltemp[NMAXLIN];  // there won't be lines longer than NMAXLIN characters
  for(;;)
  {
    (theLineNo[theCurrentFile])++;
    for(ii = 0; ii < NMAXLIN; ++ii)
    {
      ltemp[ii] = ' ';
    }
    theFiles[theCurrentFile]->getline(ltemp, NMAXLIN);

    //---------- Check for lines longer than NMAXLIN character
    for(ii = 0; ii < NMAXLIN; ++ii)
    {
      if(ltemp[ii] == '\0')
      {
        break;
      }
    }
    if(ii == NMAXLIN - 1)
    {
      ErrorInLine();
      G4String ErrMessage = "Too long line. Please split it " +
                            G4String("putting a '\\' at the end!");
      G4Exception("G4tgrFileIn::GetWordsInLine()", "InvalidInput",
                  FatalException, ErrMessage);
    }

    //---------- End of file
    if(EndOfFile())
    {
      return 0;
    }

    //---------- Convert line read to istrstream to split it in words
    std::istringstream istr_line(ltemp);

    //--------- Count how many words are there in ltemp
    //          this shouln't be needed, but SUN compiler has problems...
    G4int NoWords = 0;
    char* tt      = ltemp;

    G4String stemp(ltemp);
    do
    {
      if(*tt != ' ' && *(tt) != '\0')
      {
        if(tt == ltemp)
        {
          ++NoWords;
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 3)
          {
            G4cout << "G4tgrFileIn::GetWordsInLine() - NoWords" << NoWords
                   << ltemp << G4endl;
          }
#endif
        }
        else if(*(tt - 1) == ' ' || *(tt - 1) == '\015' || *(tt - 1) == '\t')
        {
          ++NoWords;
#ifdef G4VERBOSE
          if(G4tgrMessenger::GetVerboseLevel() >= 3)
          {
            G4cout << "G4tgrFileIn::GetWordsInLine() - NoWords" << NoWords
                   << ltemp << G4endl;
          }
#endif
        }
      }
      ++tt;
    } while((*tt != '\0') && (stemp.length() != 0));

    if(stemp.length() == 0)
    {
      NoWords = 0;
    }

    //--------- Read words from istr_line and write them into wordlist
    for(ii = 0; ii < NoWords; ++ii)
    {
      stemp = "";
      istr_line >> stemp;
      if(stemp.length() == 0)
      {
        break;
      }
      G4int comment = (G4int)stemp.find(G4String("//"));
#ifdef G4VERBOSE
      if(G4tgrMessenger::GetVerboseLevel() >= 3)
      {
        G4cout << "!!!COMMENT" << comment << stemp.c_str() << G4endl;
      }
#endif
      if(comment == 0)
      {
        break;
      }
      else if(comment > 0)
      {
        stemp = stemp.substr(0, comment);
        wordlist.push_back(stemp);
        break;
      }
      wordlist.push_back(stemp);
    }

    // These two algorithms should be the more STL-like way, but they don't
    // work for files whose lines end without '\015'=TAB (STL problem: doesn't
    // find end of string??):
    // istream_iterator<G4String, ptrdiff_t> G4String_iter(istr_line);
    // istream_iterator<G4String, ptrdiff_t> eosl;
    // copy(G4String_iter, eosl, back_inserter(wordlist));
    // typedef istream_iterator<G4String, ptrdiff_t> G4String_iter;
    // copy(G4String_iter(istr_line), G4String_iter(), back_inserter(wordlist));

    if(wordlist.size() != 0)
    {
      if((*(wordlist.end() - 1)).compare("\\") == 0)  // use '\' to mark
      {                                               // continuing line
        wordlist.pop_back();
      }
      else
      {
        break;
      }
    }
  }

  //--------- A pair of double quotes delimits a word, therefore, look for the
  //          case where there is more than one word between two double quotes
  std::vector<G4String> wordlist2;
  G4String wordq      = "";
  unsigned int imerge = 0;
  for(std::size_t jj = 0; jj < wordlist.size(); ++jj)
  {
    if(wordlist[jj].substr(0, 1) == "\"")
    {
      imerge = 1;
    }
    if(wordlist[jj][G4int(wordlist[jj].size() - 1)] == '\"')
    {
      if(imerge != 1)
      {
        G4String err1 = " word with trailing '\"' while there is no";
        G4String err2 = " previous word with leading '\"' in line ";
        G4String err  = err1 + err2;
        DumpException(err);
      }
      imerge = 2;
    }
    if(imerge == 0)
    {
      wordlist2.push_back(wordlist[jj]);
    }
    else if(imerge == 1)
    {
      if(wordq == "")
      {
        wordq.append(wordlist[jj].substr(1, wordlist[jj].size()));
      }
      else
      {
        wordq.append(wordlist[jj].substr(0, wordlist[jj].size()));
      }
      wordq.append(" ");
    }
    else if(imerge == 2)
    {
      if(wordq == "")
      {
        wordq.append(wordlist[jj].substr(1, wordlist[jj].size() - 2));
      }
      else
      {
        wordq.append(wordlist[jj].substr(0, wordlist[jj].size() - 1));
      }
      wordlist2.push_back(wordq);
      wordq  = "";
      imerge = 0;
    }
  }
  if(imerge == 1)
  {
    G4String err1 = " word with leading '\"' in line while there is no";
    G4String err2 = " later word with trailing '\"' in line ";
    G4String err  = err1 + err2;
    DumpException(err);
  }

  wordlist = wordlist2;

  // Or why not like this (?):
  // typedef std::istream_iterator<G4String, ptrdiff_t> string_iter;
  // std::copy(string_iter(istr_line), string_iter(), back_inserter(wordlist));

  // check if including a new file
  if(wordlist[0] == "#include")
  {
    if(wordlist.size() != 2)
    {
      ErrorInLine();
      G4String ErrMessage =
        "'#include' should have as second argument, the filename !";
      G4Exception("G4tgrFileIn::GetWordsInLine()", "InvalidInput",
                  FatalException, ErrMessage);
    }

#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 3)
    {
      G4cout << " G4tgrFileIn::GetWordsInLine() - Include found !" << G4endl;
    }
#endif
    OpenNewFile(wordlist[1].c_str());
    isok = GetWordsInLine(wordlist);
  }

  return isok;
}

// --------------------------------------------------------------------
void G4tgrFileIn::ErrorInLine()
{
  G4cerr << "!! EXITING: ERROR IN LINE No " << theLineNo[theCurrentFile]
         << " file: " << theNames[theCurrentFile] << " : ";
}

// --------------------------------------------------------------------
G4bool G4tgrFileIn::EndOfFile()
{
  G4bool isok = theFiles[theCurrentFile]->eof();
  if(isok)
  {
#ifdef G4VERBOSE
    if(G4tgrMessenger::GetVerboseLevel() >= 3)
    {
      G4cout << " G4tgrFileIn::EndOfFile() - EOF: " << theCurrentFile << G4endl;
    }
#endif
    --theCurrentFile;
    if(theCurrentFile != -1)  // Last file will be closed by the user
    {
      Close();
    }
  }

  // Only real closing if all files are closed
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 3)
  {
    G4cout << " G4tgrFileIn::EndOfFile() - EOF: " << isok << " "
           << theCurrentFile << G4endl;
  }
#endif
  if(theCurrentFile != -1)
  {
    return false;
  }
  else
  {
    return isok;
  }
}

// --------------------------------------------------------------------
void G4tgrFileIn::Close()
{
#ifdef G4VERBOSE
  if(G4tgrMessenger::GetVerboseLevel() >= 3)
  {
    G4cout << "G4tgrFileIn::Close() - " << theCurrentFile << ", size "
           << theFiles.size() << G4endl;
  }
#endif

  theFiles[theCurrentFile + 1]->close();
  theFiles.pop_back();
}

// --------------------------------------------------------------------
void G4tgrFileIn::DumpException(const G4String& sent)
{
  G4String Err1 = sent + " in file " + theName;
  G4String Err2 =
    " line No: " + G4UIcommand::ConvertToString(theLineNo[theCurrentFile]);
  G4String ErrMessage = Err1;
  G4Exception("G4tgrFileIn::DumpException()", "FileError", FatalException,
              ErrMessage);
}
