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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// J. Comput. Phys. 274 (2014) 841-882
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Author: Mathieu Karamitros
//
// $Id$
//
/// \file CommandLineParser.hh
/// \brief Definition of the CommandLineParser class

#ifndef COMMANDLINEPARSER_HH
#define COMMANDLINEPARSER_HH

#include "globals.hh"
#include <map>

namespace G4DNAPARSER
{
class Command
{
public:
  enum Type
  {
    WithOption,
    WithoutOption,
    OptionNotCompulsory
  };

  const G4String& GetOption() {return fOption;}
  Command::Type GetType() {return fType;}
  G4bool IsActive() {return fActive;}
  const G4String& GetDescription() {return fDescription;}
  const G4String& GetOptionName() {return fOptionName;}

private:
  friend class CommandLineParser;
  Type fType;
  G4String fOption;
  G4bool fActive;
  G4String fDescription;
  G4String fOptionName;

  Command(Type, const G4String &description = "",
          const G4String &optionName ="optionName");
  ~Command(){;}
};

class CommandLineParser
{
    static CommandLineParser* fpInstance;
    std::map<G4String, Command*> fCommandMap;
    G4bool fOptionsWereSetup;
    G4int fMaxMarkerLength;
    G4int fMaxOptionNameLength;
    G4int fVerbose;

public:
    static CommandLineParser* GetParser();
    CommandLineParser();
    ~CommandLineParser();
    static void DeleteInstance();
    int Parse(int& argc, char **argv);
    void PrintHelp();
    bool CheckIfNotHandledOptionsExists(int& argc, char** argv);
    void CorrectRemainingOptions(int& argc, char **argv);
    void AddCommand(const G4String & marker,Command::Type,
                    const G4String& description = "",
                    const G4String& optionName = "");
    Command* FindCommand(const G4String &marker);
    Command* GetCommandIfActive(const G4String &marker);
    G4bool WereOptionsSetup(){return fOptionsWereSetup;}
};
}
#endif // PARSER_HH
