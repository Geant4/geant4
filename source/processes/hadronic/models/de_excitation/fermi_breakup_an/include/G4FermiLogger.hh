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
//
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#ifndef G4FERMILOGGER_HH
#define G4FERMILOGGER_HH

#include <G4ios.hh>

#include <memory>
#include <sstream>
#include <string>
#include <string_view>

enum class G4FermiLogLevel : int
{
  TRACE = 0,
  DEBUG = 1,
  INFO = 2,
  WARN = 3,
  ERROR = 4,
  NONE = 5,
};

class G4FermiILogger
{
  public:
    virtual void Log(const std::string_view fileName, const int line,
                     const std::string_view funcName, const G4FermiLogLevel level,
                     const std::string_view msg) = 0;

    virtual ~G4FermiILogger() = default;
};

class G4FermiStreamLogger : public G4FermiILogger
{
  public:
    G4FermiStreamLogger(std::ostream& stream) : stream_(stream) {}

    void Log(const std::string_view fileName, const int line, const std::string_view funcName,
             const G4FermiLogLevel level, const std::string_view msg) override;

  private:
    std::ostream& stream_;
};

class G4FermiLogger
{
  public:
    G4FermiLogger();

    G4FermiLogger(std::shared_ptr<G4FermiILogger>&& impl,
                  const G4FermiLogLevel level = G4FermiLogLevel::TRACE)
      : impl_(std::move(impl)), level_(level)
    {}

    G4FermiLogger(const std::shared_ptr<G4FermiILogger>& impl,
                  const G4FermiLogLevel level = G4FermiLogLevel::TRACE)
      : G4FermiLogger(std::shared_ptr<G4FermiILogger>(impl), level)
    {}

    G4FermiLogger(const G4FermiLogger& other) = default;

    G4FermiLogger& operator=(const G4FermiLogger& other) = default;

    bool ShouldLog(const G4FermiLogLevel level) const
    {
      return static_cast<int>(level) >= static_cast<int>(level_);
    }

    void Log(const std::string_view fileName, const int line, const std::string_view funcName,
             const G4FermiLogLevel level, const std::string_view msg)
    {
      impl_->Log(fileName, line, funcName, level, msg);
    }

    G4FermiLogger WithLevel(const G4FermiLogLevel level) const
    {
      return G4FermiLogger(impl_, level);
    }

    operator bool() const { return impl_ != nullptr; }

    static G4FermiLogger Default() { return G4FermiLogger(GlobalImpl, GlobalLevel); }

    static inline std::shared_ptr<G4FermiILogger> GlobalImpl =
      std::make_shared<G4FermiStreamLogger>(G4cerr);
    static inline G4FermiLogLevel GlobalLevel = G4FermiLogLevel::ERROR;

  private:
    std::shared_ptr<G4FermiILogger> impl_;
    G4FermiLogLevel level_;
};

#define FERMI_LOG_MSG(logger, level, msg)                               \
  if (logger && logger.ShouldLog(level)) {                              \
    std::ostringstream sstream;                                         \
    sstream << msg;                                                     \
    logger.Log(__FILE__, __LINE__, __FUNCTION__, level, sstream.str()); \
  }

#define FERMI_LOG_TRACE(msg) FERMI_LOG_MSG(G4FermiLogger::Default(), G4FermiLogLevel::TRACE, msg)
#define FERMI_LOG_DEBUG(msg) FERMI_LOG_MSG(G4FermiLogger::Default(), G4FermiLogLevel::DEBUG, msg)
#define FERMI_LOG_INFO(msg) FERMI_LOG_MSG(G4FermiLogger::Default(), G4FermiLogLevel::INFO, msg)
#define FERMI_LOG_WARN(msg) FERMI_LOG_MSG(G4FermiLogger::Default(), G4FermiLogLevel::WARN, msg)
#define FERMI_LOG_ERROR(msg) FERMI_LOG_MSG(G4FermiLogger::Default(), G4FermiLogLevel::ERROR, msg)

#if defined(WIN32) || defined(__MINGW32__)
#  define FERMI_UNLIKELY(x) (x)
#  define FERMI_LIKELY(x) (x)
#else
#  define FERMI_UNLIKELY(x) __builtin_expect((x), 0)  // gcc/clang extension - not portable
#  define FERMI_LIKELY(x) __builtin_expect((x), 1)
#endif

#define FERMI_ASSERT_MSG(COND, MSG)                                                             \
  if (FERMI_UNLIKELY(!(COND))) {                                                                \
    std::ostringstream sstream;                                                                 \
    sstream << "assertion failed: \"" << #COND << '\"' << " at " << __FILE__ << ':' << __LINE__ \
            << '\n'                                                                             \
            << MSG;                                                                             \
    throw std::runtime_error(sstream.str());                                                    \
  }

#endif  // G4FERMILOGGER_HH
