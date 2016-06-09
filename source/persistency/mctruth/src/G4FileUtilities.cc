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
// File: G4FileUtilities.cc
//
// History:
//   01.08.24  Youhei Morita  Initial creation

#ifndef WIN32

#include "G4FileUtilities.hh"

// Implementation of Default Constructor
G4FileUtilities::G4FileUtilities()
{
}

// Implementation of Default Destructor
G4FileUtilities::~G4FileUtilities()
{
}

// Implementation of FileExists
G4bool G4FileUtilities::FileExists(const std::string file)
{
  char* c = (char *) file.c_str();

  G4int fd = ::open( c, O_RDONLY ); 
  // G4int error = errno;
  if ( fd != -1 ) {
    ::close( fd );
    return true;
  } else {
    return false;
  }
}

// Implementation of CopyFile
G4int G4FileUtilities::CopyFile(const std::string srcFile, const std::string dstFile)
{
  std::string cmd = "cp " + srcFile + " " + dstFile;
  return Shell( cmd );
}

// Implementation of DeleteFile
G4int G4FileUtilities::DeleteFile(const std::string file, const std::string option)
{
  std::string cmd = "rm " + option + " " + file;
  return Shell( cmd );
}

#endif
