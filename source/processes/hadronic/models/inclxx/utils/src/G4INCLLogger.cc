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
// INCL++ intra-nuclear cascade model
// Alain Boudard, CEA-Saclay, France
// Joseph Cugnon, University of Liege, Belgium
// Jean-Christophe David, CEA-Saclay, France
// Pekka Kaitaniemi, CEA-Saclay, France, and Helsinki Institute of Physics, Finland
// Sylvie Leray, CEA-Saclay, France
// Davide Mancusi, CEA-Saclay, France
//
#define INCLXX_IN_GEANT4_MODE 1

#include "globals.hh"

#include "G4INCLLogger.hh"
#include "G4INCLGlobals.hh"

#ifdef INCLXX_IN_GEANT4_MODE
#include <cstdlib>
#endif

namespace G4INCL {
#if defined(INCL_DEBUG_LOG) && !defined(INCLXX_IN_GEANT4_MODE)
  std::string typeToString(const MessageType t) {
    if(t == ErrorMsg)
      return std::string("Error");
    else if(t == FatalMsg)
      return std::string ("Fatal");
    else if(t == WarningMsg)
      return std::string("Warning");
    else if(t == DebugMsg)
      return std::string("Debug");
    else if(t == InfoMsg)
      return std::string("");
    else if(t == DataBlockMsg)
      return std::string("DataBlock");
    else
      return std::string("Unknown");
  }

  void LoggerSlave::logMessage(const MessageType type, const std::string &fileName, const G4int lineNumber, std::string const &s, const G4bool prefixHash) const {
    std::string cont = (prefixHash ? "\n# " : "\n");
    std::string header = (prefixHash ? "# " : "");
    if(type!=InfoMsg) {
      std::stringstream headerss;
      headerss << typeToString(type) << " [";
      cont += headerss.str();
      headerss <<
        fileName.substr(fileName.find_last_of("/")+1) <<
        ":" << lineNumber << "] ";
      header += headerss.str();
      cont.append(header.size() - cont.size() - 1, '.');
      cont += "] ";
    }

    std::string message(s);
    String::replaceAll(message, "\n", cont, s.size()-2);
    (*logStream) << header << message;

    return;
  }

  void LoggerSlave::logDataBlock(const std::string &block, const std::string &fileName, const G4int lineNumber) const {
    (*logStream) << typeToString(DataBlockMsg) << " [" <<
      fileName.substr(fileName.find_last_of("/")+1) <<
      ":" << lineNumber << "] " << '\n'
      << "BEGINDATA"
      << block
      << "ENDDATA" << '\n';
  }

  namespace Logger {

    namespace {
      G4ThreadLocal LoggerSlave *theLoggerSlave = NULL;
    }

    void logMessage(const MessageType type, std::string const &fileName, const G4int lineNumber, std::string const &s, const G4bool prefixHash) {
      if(theLoggerSlave)
        theLoggerSlave->logMessage(type, fileName, lineNumber, s, prefixHash);
    }

    void flush() {
      if(theLoggerSlave)
        theLoggerSlave->flush();
    }

    void dataBlock(const std::string &block, const std::string &fileName, const G4int lineNumber) {
      theLoggerSlave->logDataBlock(block, fileName, lineNumber);
    }

    void setLoggerSlave(LoggerSlave * const logger) { theLoggerSlave = logger; }

    void setVerbosityLevel(G4int lvl) {
      if(theLoggerSlave)
        theLoggerSlave->setVerbosityLevel(lvl);
    }

    G4int getVerbosityLevel() {
      if(theLoggerSlave)
        return theLoggerSlave->getVerbosityLevel();
      else
        return 0;
    }

    void deleteLoggerSlave() {
      delete theLoggerSlave;
      theLoggerSlave=NULL;
    }

    void initialize(Config const * const theConfig) {
      setLoggerSlave(new LoggerSlave(theConfig->getLogFileName(), theConfig->getVerbosity()));
    }

  }

#else // defined(INCL_DEBUG_LOG) && !defined(INCLXX_IN_GEANT4_MODE)

  namespace Logger {

    namespace {
      G4ThreadLocal G4int verbosityLevel = 0;
    }

    void initVerbosityLevelFromEnvvar() {
      const char * const envVar = getenv("G4INCL_DEBUG_VERBOSITY");
      if(envVar) {
        std::stringstream verbss(envVar);
        verbss >> verbosityLevel;
      } else {
        verbosityLevel = 0;
      }
    }

    G4int getVerbosityLevel() {
      return verbosityLevel;
    }

  }

#endif // defined(INCL_DEBUG_LOG) && !defined(INCLXX_IN_GEANT4_MODE)
}
