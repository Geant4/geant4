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

#ifndef G4INCLLogger_hh
#define G4INCLLogger_hh 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>

#ifdef INCLXX_IN_GEANT4_MODE
#include "G4ios.hh"
#endif

#include "G4INCLRandom.hh"
#include "G4INCLConfig.hh"

namespace G4INCL {

  /**
   * Verbosity scale from 0 (fatal errors only) to 10 (print everything)
   */
  enum MessageType { InfoMsg = 1,
    FatalMsg = 2,
    ErrorMsg = 3,
    WarningMsg = 4,
    DebugMsg = 7,
    DataBlockMsg = 10,
    ZeroMsg = 0 };

#if defined(INCL_DEBUG_LOG) && !defined(INCLXX_IN_GEANT4_MODE)

  class LoggerSlave {
    public:
      // By default, log fatal errors, errors and warnings
      LoggerSlave(std::string const &logFileName, const G4int verbosity=4) :
        logStream(0),
        verbosityLevel(verbosity)
    {
      if(logFileName=="-") {
        logStream = &(std::cout);
        logToStdout = true;
      } else {
        logToStdout = false;
        logStream = new std::ofstream(logFileName.c_str());
        if(!logStream)
        {
          std::cerr << "Fatal error: couldn't open log file " << logFileName << std::endl;
          std::exit(EXIT_FAILURE);
        }
      }

      // Spell out "true" and "false" when logging G4bool variables
      std::boolalpha(*logStream);
    };
      ~LoggerSlave() {
        if(!logToStdout)
          delete logStream;
      };

      /**
       * Set the verbosity level
       */
      void setVerbosityLevel(G4int lvl) { verbosityLevel = lvl; }

      /**
       * Get the verbosity level
       */
      G4int getVerbosityLevel() { return verbosityLevel; }

      /// \brief Write the log message.
      void logMessage(const MessageType type, const std::string &fileName, const G4int lineNumber, std::string const &s, const G4bool prefixHash=true) const;

      /// \brief Flush the log stream
      void flush() { logStream->flush(); }

      /// \brief Log a data block.
      void logDataBlock(const std::string &block, const std::string &fileName, const G4int lineNumber) const;

      typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
      typedef CoutType& (*StandardEndLine)(CoutType&);
      /// \brief Overload << operator to support std::endl.
      LoggerSlave const &operator<<(StandardEndLine const &manip) const {
        manip(*logStream);
        return *this;
      }

      /// \brief Overloaded << operator to provide a stream-like API.
      template<typename T>
        LoggerSlave const &operator<<(const T &t) const {
          (*logStream) << t;
          return *this;
        }

    private:
      std::ostream *logStream;
      G4int verbosityLevel;
      G4bool logToStdout;
  };

  namespace Logger {
      /// \brief Log a message.
      void logMessage(const MessageType type, std::string const &fileName, const G4int lineNumber, std::string const &s, const G4bool prefixHash=true);

      /// \brief Flush the log stream
      void flush();

      /// \brief Log a data block.
      void dataBlock(const std::string &block, const std::string &fileName, const G4int lineNumber);

      /// \brief Set the slave Logger.
      void setLoggerSlave(LoggerSlave * const logger);

      /// \brief Set the verbosity of the slave Logger.
      void setVerbosityLevel(G4int lvl);

      /// \brief Get the verbosity of the slave Logger.
      G4int getVerbosityLevel();

      /// \brief Delete the slave Logger.
      void deleteLoggerSlave();

      /// \brief Initialize the Logger.
      void initialize(Config const * const theConfig);

  }

  // Macro definitions for line numbering in log files!
#define INCL_FATAL(x) \
  if(true) {\
    std::stringstream ss_;\
    ss_ << x;\
    ss_ << "Random seeds at the beginning of this event: " << G4INCL::Random::getSavedSeeds() << std::endl;\
    G4INCL::Logger::logMessage(G4INCL::FatalMsg, __FILE__,__LINE__, ss_.str());\
    G4INCL::Logger::flush();\
    std::exit(EXIT_FAILURE);\
  } else (void)0
#define INCL_ERROR(x) \
  if(G4INCL::ErrorMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::stringstream ss_;\
    ss_ << x;\
    ss_ << "Random seeds at the beginning of this event: " << G4INCL::Random::getSavedSeeds() << std::endl;\
    G4INCL::Logger::logMessage(G4INCL::ErrorMsg, __FILE__,__LINE__, ss_.str());\
  } else (void)0
#define INCL_WARN(x) \
  if(G4INCL::WarningMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::stringstream ss_;\
    ss_ << x;\
    G4INCL::Logger::logMessage(G4INCL::WarningMsg, __FILE__,__LINE__, ss_.str());\
  } else (void)0
#define INCL_INFO(x) \
  if(G4INCL::InfoMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::stringstream ss_;\
    ss_ << x;\
    G4INCL::Logger::logMessage(G4INCL::InfoMsg, __FILE__,__LINE__, ss_.str());\
  } else (void)0
#define INCL_INFO_NOCOMMENT(x) \
  if(G4INCL::InfoMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::stringstream ss_;\
    ss_ << x;\
    G4INCL::Logger::logMessage(G4INCL::InfoMsg, __FILE__,__LINE__, ss_.str(), false);\
  } else (void)0
#define INCL_DEBUG(x) \
  if(G4INCL::DebugMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::stringstream ss_;\
    ss_ << x;\
    G4INCL::Logger::logMessage(G4INCL::DebugMsg, __FILE__,__LINE__, ss_.str());\
  } else (void)0
#define INCL_DATABLOCK(x) \
  if(G4INCL::DataBlockMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    G4INCL::Logger::dataBlock(x,__FILE__,__LINE__);\
  } else (void)0

#else // defined(INCL_DEBUG_LOG) && !defined(INCLXX_IN_GEANT4_MODE)
  namespace Logger {
    void initVerbosityLevelFromEnvvar();
    G4int getVerbosityLevel();
  }

#define INCL_FATAL(x) \
  if(true) {\
    std::stringstream ss_;\
    ss_ << x;\
    std::stringstream location_;\
    std::string fileName_(__FILE__);\
    location_ << fileName_.substr(fileName_.find_last_of("/")+1) << ":" << __LINE__;\
    G4Exception(location_.str().c_str(), "INCLXX0000", EventMustBeAborted, ss_.str().c_str());\
  } else (void)0
#define INCL_ERROR(x) \
  if(G4INCL::ErrorMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::string fileName_(__FILE__);\
    std::stringstream ss_;\
    ss_ << "INCL++ error [" << fileName_.substr(fileName_.find_last_of("/")+1) << ":" << __LINE__ << "] " << x;\
    G4cout << ss_.str() << '\n';\
  } else (void)0
#define INCL_WARN(x) \
  if(G4INCL::WarningMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::string fileName_(__FILE__);\
    std::stringstream ss_;\
    ss_ << "INCL++ warning [" << fileName_.substr(fileName_.find_last_of("/")+1) << ":" << __LINE__ << "] " << x;\
    G4cout << ss_.str() << '\n';\
  } else (void)0
#define INCL_INFO(x);
#define INCL_DEBUG(x) \
  if(G4INCL::DebugMsg <= G4INCL::Logger::getVerbosityLevel()) {\
    std::string fileName_(__FILE__);\
    std::stringstream ss_;\
    ss_ << "INCL++ debug [" << fileName_.substr(fileName_.find_last_of("/")+1) << ":" << __LINE__ << "] " << x;\
    G4cout << ss_.str() << '\n';\
  } else (void)0
#define INCL_DATABLOCK(x);

#endif // defined(INCL_DEBUG_LOG) && !defined(INCLXX_IN_GEANT4_MODE)
}
#endif
