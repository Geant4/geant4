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
// File: G4VTransactionManager.hh
//
// History:
//   01.07.18  Youhei Morita  Initial creation (with "fadsclass")

#ifndef V_TRANSACTION_MANAGER_HH
#define V_TRANSACTION_MANAGER_HH 1

#include "G4Types.hh"
#include <string>

// Class Description:
//   Abstract base class for starting and commiting or aborting data transaction

class G4VTransactionManager
{
    public: // With description
      G4VTransactionManager() : m_verbose(0) {};
      // Constructor

      virtual ~G4VTransactionManager() {};
      // Destructor

    public: // With description
      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      virtual G4bool SelectReadFile(std::string obj, std::string file)=0;
      // Set the input file name and open it for the object type "obj".

      virtual G4bool SelectWriteFile(std::string obj, std::string file)=0;
      // Set the output file name and open it for the object type "obj".

      virtual G4bool StartUpdate()=0;
      // Start an update transaction for event store and retrieve

      virtual G4bool StartRead()=0;
      // Start a read-only transaction for event store and retrieve

      virtual void Commit()=0;
      // commit the transaction for event store and retrieve

      virtual void Abort()=0;
      // abort the transaction for event store and retrieve

    protected:
      G4int m_verbose;

}; // End of class G4VTransactionManager

#endif

