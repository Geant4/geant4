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
// File: G4FileUtilities.hh
//
// History:
//   01.08.24  Youhei Morita  Initial creation

#ifndef FILE_UTILITIES_HH
#define FILE_UTILITIES_HH 1

#include "G4Types.hh"

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include "g4std/iostream"


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
      G4bool FileExists(const G4std::string file);
      // checks if the "file" exists.  returns true if it does.

      G4std::string StrErrNo() const { return ::strerror(errno); };
      // returns the error message of the last system call as string.

      int Shell(G4std::string s) { return ::system(s.c_str()); };
      // execute the shell command.  returns zero if success.

      int CopyFile(const G4std::string srcFile, const G4std::string dstFile);
      // copies the "srcFile" to "dstFile".  returns zero if success.

      int DeleteFile(const G4std::string file, const G4std::string option);
      // deletes the "file" with the "option".  returns zero if success.

      G4std::string GetEnv(const G4std::string env) { return ::getenv(env.c_str()); };
      // retuns the value of environment variable as string.

}; // End of class G4FileUtilities

#endif

