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
// G4FileUtilities implementation
//
// Author: Youhei Morita, 24.08.2001
// --------------------------------------------------------------------

#ifndef WIN32

#include "G4FileUtilities.hh"

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

// --------------------------------------------------------------------
G4FileUtilities::G4FileUtilities()
{
}

// --------------------------------------------------------------------
G4FileUtilities::~G4FileUtilities()
{
}

// --------------------------------------------------------------------
G4bool G4FileUtilities::FileExists(const G4String& file)
{
  char* c = (char*) file.c_str();

  G4int fd = ::open(c, O_RDONLY);
  // G4int error = errno;
  if(fd != -1)
  {
    ::close(fd);
    return true;
  }
  else
  {
    return false;
  }
}

// --------------------------------------------------------------------
G4String G4FileUtilities::StrErrNo() const
{
  return ::strerror(errno);
}

// --------------------------------------------------------------------
G4int G4FileUtilities::Shell(const G4String& str)
{
  return ::system(str.c_str());
}

// --------------------------------------------------------------------
G4int G4FileUtilities::CopyFile(const G4String& srcFile,
                                const G4String& dstFile)
{
  G4String cmd = "cp " + srcFile + " " + dstFile;
  return Shell(cmd);
}

// --------------------------------------------------------------------
G4int G4FileUtilities::DeleteFile(const G4String& file,
                                  const G4String& option)
{
  G4String cmd = "rm " + option + " " + file;
  return Shell(cmd);
}

// --------------------------------------------------------------------
G4String G4FileUtilities::GetEnv(const G4String& env)
{
  return ::getenv(env.c_str());
}

#endif
