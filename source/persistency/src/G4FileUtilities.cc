// $Id: G4FileUtilities.cc,v 1.2 2002-12-04 10:25:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4FileUtilities.cc
//
// History:
//   01.08.24  Youhei Morita  Initial creation

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

// End of G4FileUtilities.cc

