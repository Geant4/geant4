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
#include <iostream>


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
      G4bool FileExists(const std::string file);
      // checks if the "file" exists.  returns true if it does.

      std::string StrErrNo() const { return ::strerror(errno); };
      // returns the error message of the last system call as string.

      G4int Shell(std::string str) { return ::system(str.c_str()); };
      // execute the shell command.  returns zero if success.

      G4int CopyFile(const std::string srcFile, const std::string dstFile);
      // copies the "srcFile" to "dstFile".  returns zero if success.

      G4int DeleteFile(const std::string file, const std::string option);
      // deletes the "file" with the "option".  returns zero if success.

      std::string GetEnv(const std::string env) { return ::getenv(env.c_str()); };
      // retuns the value of environment variable as string.

}; // End of class G4FileUtilities

#endif

