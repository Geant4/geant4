// $Id: G4FileUtilities.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4FileUtilities.hh
//
// History:
//   01.08.24  Youhei Morita  Initial creation

#ifndef FILE_UTILITIES_HH
#define FILE_UTILITIES_HH 1

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <iostream.h>

// Class Description:
//   File utilities to access files with POSIX interface.

class G4FileUtilities
{
    public: // With description
      G4FileUtilities();
      // Constructor

      ~G4FileUtilities();
      // Destructor

    public: // With description
      bool FileExists(const std::string file);
      // checks if the "file" exists.  returns true if it does.

      std::string StrErrNo() const { return ::strerror(errno); };
      // returns the error message of the last system call as string.

      int Shell(std::string s) { return ::system(s.c_str()); };
      // execute the shell command.  returns zero if success.

      int CopyFile(const std::string srcFile, const std::string dstFile);
      // copies the "srcFile" to "dstFile".  returns zero if success.

      int DeleteFile(const std::string file, const std::string option);
      // deletes the "file" with the "option".  returns zero if success.

      std::string GetEnv(const std::string env) { return ::getenv(env.c_str()); };
      // retuns the value of environment variable as string.

}; // End of class G4FileUtilities

#endif

