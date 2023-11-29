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
// Global environment utility functions:
//
// G4GetEnv<T>
//      Simplifies getting environment variables
//      Automatic conversion to non-string types
//      Records the values used from the environment
// G4GetDataEnv
//      For data library paths
//      Will issue a G4Exception if not set
// G4PrintEnv
//      Provide a way for users to determine (and log) the environment
//      variables were used as settings in simulation

// Author: Jonathan Madsen, 25 October 2018
// ---------------------------------------------------------------------------
#ifndef G4ENVIRONMENTUTILS_HH
#define G4ENVIRONMENTUTILS_HH

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>

#include "G4Exception.hh"
#include "G4ExceptionSeverity.hh"
#include "G4String.hh"
#include "G4ios.hh"

class G4EnvSettings
{
  // Static singleton class storing environment variables and
  // their values that were used by Geant4 in the simulation

 public:
  using string_t   = std::string;
  using env_map_t  = std::map<string_t, string_t>;
  using env_pair_t = std::pair<string_t, string_t>;

  static G4EnvSettings* GetInstance()
  {
    static auto* _instance = new G4EnvSettings();
    return _instance;
  }

  template <typename _Tp>
  void insert(const std::string& env_id, _Tp val)
  {
    std::stringstream ss;
    ss << val;
    // lock for MT mode, use C++ type not Geant4 because this file
    // is included by the those headers
    static std::mutex _mutex;
    _mutex.lock();
    m_env.insert(env_pair_t(env_id, ss.str()));
    _mutex.unlock();
  }

  const env_map_t& get() const { return m_env; }

  friend std::ostream& operator<<(std::ostream& os, const G4EnvSettings& env)
  {
    std::stringstream filler;
    filler.fill('#');
    filler << std::setw(90) << "";
    std::stringstream ss;
    ss << filler.str() << "\n# Environment settings:\n";
    for(const auto& itr : env.get())
    {
      ss << "# " << std::setw(35) << std::right << itr.first << "\t = \t"
         << std::left << itr.second << "\n";
    }
    ss << filler.str();
    os << ss.str() << std::endl;
    return os;
  }

 private:
  env_map_t m_env;
};

// ---------------------------------------------------------------------------
//  Use this function to get an environment variable setting +
//  a default if not defined, e.g.
//      int num_threads =
//          G4GetEnv<int>("G4FORCENUMBEROFTHREADS",
//                        std::thread::hardware_concurrency());
template <typename _Tp>
_Tp G4GetEnv(const std::string& env_id, _Tp _default = _Tp())
{
  char* env_var = std::getenv(env_id.c_str());
  if(env_var)
  {
    std::string str_var = std::string(env_var);
    std::istringstream iss(str_var);
    _Tp var = _Tp();
    iss >> var;
    // record value defined by environment
    G4EnvSettings::GetInstance()->insert<_Tp>(env_id, var);
    return var;
  }
  // record default value
  G4EnvSettings::GetInstance()->insert<_Tp>(env_id, _default);

  // return default if not specified in environment
  return _default;
}

// ---------------------------------------------------------------------------
//  Use this function to get an environment variable setting +
//  a default if not defined, e.g.
//      int num_threads =
//          GetEnv<int>("FORCENUMBEROFTHREADS",
//                      std::thread::hardware_concurrency());
template <>
inline G4bool G4GetEnv(const std::string& env_id, bool _default)
{
  char* env_var = std::getenv(env_id.c_str());
  if(env_var != nullptr)
  {
    // record value defined by environment
    G4EnvSettings::GetInstance()->insert<bool>(env_id, true);
    return true;
  }
  // record default value
  G4EnvSettings::GetInstance()->insert<bool>(env_id, false);

  // return default if not specified in environment
  return _default;
}

// ---------------------------------------------------------------------------
//  Use this function to get an environment variable setting +
//  a default if not defined and a message about the setting, e.g.
//      int num_threads =
//          G4GetEnv<int>("G4FORCENUMBEROFTHREADS",
//                        std::thread::hardware_concurrency(),
//                        "Forcing number of threads");
template <typename _Tp>
_Tp G4GetEnv(const std::string& env_id, _Tp _default, const std::string& msg)
{
  char* env_var = std::getenv(env_id.c_str());
  if(env_var)
  {
    std::string str_var = std::string(env_var);
    std::istringstream iss(str_var);
    _Tp var = _Tp();
    iss >> var;
    G4cout << "Environment variable \"" << env_id << "\" enabled with "
           << "value == " << var << ". " << msg << G4endl;
    // record value defined by environment
    G4EnvSettings::GetInstance()->insert<_Tp>(env_id, var);
    return var;
  }
  // record default value
  G4EnvSettings::GetInstance()->insert<_Tp>(env_id, _default);

  // return default if not specified in environment
  return _default;
}

// ---------------------------------------------------------------------------
//  Use this function to get a data directory environment variable setting +
//  and raise a G4Exception if the value is not set, e.g.
//
//      G4String filename = G4GetDataEnv("G4ENSDFSTATEDATA",
//                                       "G4NuclideTable", "PART70000",
//                                       FatalException,
//                                       "G4ENSDFSTATEDATA environment variable"
//                                       " must be set");
inline G4String G4GetDataEnv(const std::string& env_id,
                             const char* originOfException,
                             const char* exceptionCode,
                             G4ExceptionSeverity severity,
                             const char* description)
{
  char* env_var = std::getenv(env_id.c_str());
  if(env_var != nullptr)
  {
    std::string str_var = std::string(env_var);
    std::istringstream iss(str_var);
    G4String var = "";
    iss >> var;
    // record value defined by environment
    G4EnvSettings::GetInstance()->insert<G4String>(env_id, var);
    return var;
  }

  // issue an exception
  G4Exception(originOfException, exceptionCode, severity, description);

  // return default initialized
  return "";
}

const char* G4FindDataDir(const char*);

// ---------------------------------------------------------------------------
// Use this function to print the environment
//
inline void G4PrintEnv(std::ostream& os = G4cout)
{
  os << (*G4EnvSettings::GetInstance());
}

#endif /* G4ENVIRONMENTUTILS_HH */
