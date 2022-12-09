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
// G4UIbatch
//
// Author: M.Asai, 2000
// --------------------------------------------------------------------

#include "G4UIbatch.hh"
#include "G4UImanager.hh"
#include <vector>
#include <string>

// --------------------------------------------------------------------
static void Tokenize(const G4String& str, std::vector<G4String>& tokens)
{
  const char* delimiter = " ";

  G4String::size_type pos0 = str.find_first_not_of(delimiter);
  G4String::size_type pos  = str.find_first_of(delimiter, pos0);

  while(pos != G4String::npos || pos0 != G4String::npos)
  {
    if(str[(G4int)pos0] == '\"')
    {
      pos = str.find_first_of("\"", pos0 + 1);
      if(pos != G4String::npos)
      {
        pos++;
      }
    }
    if(str[(G4int)pos0] == '\'')
    {
      pos = str.find_first_of("\'", pos0 + 1);
      if(pos != G4String::npos)
      {
        pos++;
      }
    }

    tokens.emplace_back(str.substr(pos0, pos - pos0));
    pos0 = str.find_first_not_of(delimiter, pos);
    pos  = str.find_first_of(delimiter, pos0);
  }
}

// --------------------------------------------------------------------
G4UIbatch::G4UIbatch(const char* fileName, G4UIsession* prevSession)
  : G4UIsession(1)
  , previousSession(prevSession)
{
  macroStream.open(fileName, std::ios::in);
  if(macroStream.fail())
  {
    G4cerr << "ERROR: Can not open a macro file <" << fileName
           << ">. Set macro path with \"/control/macroPath\" if needed."
           << G4endl;
    lastRC = fParameterUnreadable;
  }
  else
  {
    isOpened = true;
  }

  G4UImanager::GetUIpointer()->SetSession(this);
}

// --------------------------------------------------------------------
G4UIbatch::~G4UIbatch()
{
  if(isOpened)
  {
    macroStream.close();
  }
}

// --------------------------------------------------------------------
G4String G4UIbatch::ReadCommand()
{
  enum
  {
    BUFSIZE = 4096
  };
  static G4ThreadLocal char* linebuf = nullptr;
  if(linebuf == nullptr)
  {
    linebuf = new char[BUFSIZE];
  }
  const char ctrM = 0x0d;

  G4String cmdtotal = "";
  G4bool qcontinued = false;
  while(macroStream.good())
  {
    macroStream.getline(linebuf, BUFSIZE);

    G4String cmdline(linebuf);

    // TAB-> ' ' conversion
    G4String::size_type nb = 0;
    while((nb = cmdline.find('\t', nb)) != G4String::npos)
    {
      cmdline.replace(nb, 1, " ");
    }

    // strip
    G4StrUtil::strip(cmdline);
    G4StrUtil::rstrip(cmdline, ctrM);

    // skip null line if single line
    if(!qcontinued && cmdline.empty())
    {
      continue;
    }

    // '#' is treated as echoing something
    if(cmdline[(std::size_t) 0] == '#')
    {
      return cmdline;
    }

    // tokenize...
    std::vector<G4String> tokens;
    Tokenize(cmdline, tokens);
    qcontinued = false;
    for(G4int i = 0; i < G4int(tokens.size()); ++i)
    {
      // string after '#" is ignored
      if(tokens[i][(std::size_t) 0] == '#')
      {
        break;
      }
      // '\' or '_' is treated as continued line.
      if(tokens[i] == "\\" || tokens[i] == "_")
      {
        qcontinued = true;
        // check nothing after line continuation character
        if(i != G4int(tokens.size()) - 1)
        {
          G4Exception("G4UIbatch::ReadCommand", "UI0003", JustWarning,
                      "unexpected character after line continuation character");
        }
        break;  // stop parsing
      }
      cmdtotal += tokens[i];
      cmdtotal += " ";
    }

    if(qcontinued)
    {
      continue;  // read the next line
    }

    if(!cmdtotal.empty())
    {
      break;
    }
    if(macroStream.eof())
    {
      break;
    }
  }

  // strip again
  G4StrUtil::strip(cmdtotal);

  // finally,
  if(macroStream.eof() && cmdtotal.empty())
  {
    return "exit";
  }

  return cmdtotal;
}

// --------------------------------------------------------------------
G4int G4UIbatch::ExecCommand(const G4String& command)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  G4int rc        = UI->ApplyCommand(command);

  switch(rc)
  {
    case fCommandSucceeded:
      break;
    case fCommandNotFound:
      G4cerr << "***** COMMAND NOT FOUND <" << command << "> *****" << G4endl;
      break;
    case fIllegalApplicationState:
      G4cerr << "***** Illegal application state <" << command << "> *****"
             << G4endl;
      break;
    default:
      G4int pn = rc % 100;
      G4cerr << "***** Illegal parameter (" << pn << ") <" << command
             << "> *****" << G4endl;
  }

  return rc;
}

// --------------------------------------------------------------------
G4UIsession* G4UIbatch::SessionStart()
{
  if(!isOpened)
  {
    return previousSession;
  }

  while(true)
  {
    G4String newCommand = ReadCommand();

    if(newCommand == "exit")
    {
      break;
    }

    // just echo something
    if(newCommand[(std::size_t) 0] == '#')
    {
      if(G4UImanager::GetUIpointer()->GetVerboseLevel() == 2)
      {
        G4cout << newCommand << G4endl;
      }
      continue;
    }

    // execute command
    G4int rc = ExecCommand(newCommand);
    if(rc != fCommandSucceeded)
    {
      G4cerr << G4endl << "***** Batch is interrupted!! *****" << G4endl;
      lastRC = rc;
      break;
    }
  }

  return previousSession;
}

// --------------------------------------------------------------------
void G4UIbatch::PauseSessionStart(const G4String& Prompt)
{
  G4cout << "Pause session <" << Prompt << "> start." << G4endl;

  SessionStart();

  G4cout << "Pause session <" << Prompt << "> Terminate." << G4endl;
}
