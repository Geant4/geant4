// $Id: G4VDCIOentry.hh,v 1.1 2002-11-24 13:45:23 morita Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// File: G4VDCIOentry.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef VDCIO_ENTRY_T_HH
#define VDCIO_ENTRY_T_HH 1

#include <string>
#include "G4PersistencyCenter.hh"

// Class Description:
//   Abstract base class for digits collection I/O manager entry

class G4VDCIOentry
{
    public: // With description
      G4VDCIOentry(std::string n);
      // Constructor

      virtual ~G4VDCIOentry() {};
      // Destructor

    public: // With description
      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      std::string GetName() { return m_name; };
      // Returns the name of the DC I/O manager entry

      virtual void CreateDCIOmanager(std::string detName, std::string colName) {};
      // virtual method for creating DC I/O manager for the detector

    protected:
      int m_verbose;

    private:
      std::string m_name;

}; // End of class G4VDCIOentry

#endif

