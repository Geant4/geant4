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
// File: G4VHCIOentry.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef VHCIO_ENTRY_T_HH
#define VHCIO_ENTRY_T_HH 1

#include "G4Types.hh"
#include <string>
#include "G4PersistencyCenter.hh"

// Class Description:
//   Abstract base class for hits collection I/O manager entry

class G4VHCIOentry
{
    public: // With description
      G4VHCIOentry(G4std::string n);
      // Constructor

      virtual ~G4VHCIOentry() {};
      // Destructor

    public: // With description
      void SetVerboseLevel(int v) { m_verbose = v; };
      // Set verbose level.

      G4std::string GetName() { return m_name; };
      // Returns the name of the HC I/O manager entry

      virtual void CreateHCIOmanager(G4std::string detName, G4std::string colName) {};
      // virtual method for creating HC I/O manager for the detector

    protected:
      int m_verbose;

    private:
      G4std::string m_name;

}; // End of class G4VHCIOentry

#endif

