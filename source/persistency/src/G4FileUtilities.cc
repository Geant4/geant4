// $Id: G4FileUtilities.cc,v 1.1 2002-11-24 13:45:24 morita Exp $
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
bool G4FileUtilities::FileExists(const std::string file)
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
int G4FileUtilities::CopyFile(const std::string srcFile, const std::string dstFile)
{
  std::string cmd = "cp " + srcFile + " " + dstFile;
  return Shell( cmd );
}

// Implementation of DeleteFile
int G4FileUtilities::DeleteFile(const std::string file, const std::string option)
{
  std::string cmd = "rm " + option + " " + file;
  return Shell( cmd );
}

// End of G4FileUtilities.cc

