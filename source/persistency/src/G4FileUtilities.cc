//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
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
G4bool G4FileUtilities::FileExists(const G4std::string file)
{
  char* c = (char *) file.c_str();

  int fd = ::open( c, O_RDONLY ); 
  // int error = errno;
  if ( fd != -1 ) {
    ::close( fd );
    return true;
  } else {
    return false;
  }
}

// Implementation of CopyFile
int G4FileUtilities::CopyFile(const G4std::string srcFile, const G4std::string dstFile)
{
  G4std::string cmd = "cp " + srcFile + " " + dstFile;
  return Shell( cmd );
}

// Implementation of DeleteFile
int G4FileUtilities::DeleteFile(const G4std::string file, const G4std::string option)
{
  G4std::string cmd = "rm " + option + " " + file;
  return Shell( cmd );
}

#endif
