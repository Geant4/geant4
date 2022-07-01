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

#include <mutex>

#ifndef G4GMAKE
#include "G4FindDataDir.hh"
#include "G4Filesystem.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>

#if defined(_MSC_VER)
#define setenv(name, value, overwrite) _putenv_s(name, value)
#endif

using namespace G4fs;

static const char * const system_paths[] = {
  GEANT4_INSTALL_FULL_DATADIR,
  CMAKE_INSTALL_PREFIX,
#if defined(_MSC_VER)
  "C:\\Program Files",
  "C:\\Geant4"
#else
  "/usr/local",
  "/usr",
  "/cvmfs/geant4.cern.ch"
#endif
};

static const char * const data_paths[] = {
  ".",
  GEANT4_INSTALL_DATADIR,
  CMAKE_INSTALL_DATADIR,
  "share/Geant4/data",
  "share/geant4/data",
  "share/data",
  "data"
};

static const char* G4GetDataDir(const char* name)
{
  for (const auto& d : dataset_definitions)
    if (strcmp(name, d.env) == 0)
      return d.dir;

  return nullptr;
}

static const char* G4FindDataDir(const char* name, const path& prefix, const path& dataset)
{
  if (!is_directory(prefix))
    return nullptr;

  for (const auto data_path : data_paths) {
    path datadir = prefix;
    if (strcmp(data_path,".") == 0)
      datadir /=  dataset;
    else
      datadir /=  path(data_path) / dataset;
    if (is_directory(absolute(datadir)))
      return setenv(name, absolute(datadir).string().c_str(), 0) == 0 ? std::getenv(name) : nullptr;
  }

  return nullptr;
}

const char* G4FindDataDir(const char* name)
{ 
#if defined(G4MULTITHREADED)
  static std::mutex mutex;
  std::lock_guard<std::mutex> lock(mutex);
#endif

  /* If environment variable is set for this dataset, use it */
  if (const char *datadir = std::getenv(name))
    return datadir;

  /* If we know which directory/version to search for, try to find it */
  if (const char *dataset = G4GetDataDir(name)) {
    /* If GEANT4_DATA_DIR environment variable is set, use it and don't search further */
    if (const char *basedir = std::getenv("GEANT4_DATA_DIR"))
      return G4FindDataDir(name, basedir, dataset);

    /* If GEANT4_DATA_DIR environment variable is not set, search in default system paths */
    for (const auto prefix : system_paths)
      if (const auto datadir = G4FindDataDir(name, prefix, dataset))
        return datadir;
  }

  return nullptr;
}

#else
// Support GNUmake builds
#include <cstdlib>

const char* G4FindDataDir(const char* name)
{
#if defined(G4MULTITHREADED)
  static std::mutex mutex;
  std::lock_guard<std::mutex> lock(mutex);
#endif
  /* If environment variable is set for this dataset, use it */
  if (const char *datadir = std::getenv(name))
    return datadir;

  return nullptr;
}
#endif
