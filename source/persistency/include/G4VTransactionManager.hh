// $Id: G4VTransactionManager.hh,v 1.2 2002-12-04 10:25:49 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

      virtual G4bool SelectReadFile(G4std::string obj, G4std::string file)=0;
      // Set the input file name and open it for the object type "obj".

      virtual G4bool SelectWriteFile(G4std::string obj, G4std::string file)=0;
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

