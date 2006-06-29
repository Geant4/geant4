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
// File: G4VDCIOentry.hh
//
// History:
//   '01.09.12  Youhei Morita  Initial creation

#ifndef VDCIO_ENTRY_T_HH
#define VDCIO_ENTRY_T_HH 1

#include "G4Types.hh"
#include <string>
#include "G4PersistencyCenter.hh"

// Class Description:
//   Abstract base class for digits collection I/O manager entry

class G4VDCIOentry
{
    public: // With description
      G4VDCIOentry(std::string n);
      // Constructor

      virtual ~G4VDCIOentry() {}
      // Destructor

    public: // With description
      void SetVerboseLevel(G4int v) { m_verbose = v; }
      // Set verbose level.

      std::string GetName() { return m_name; }
      // Returns the name of the DC I/O manager entry

      virtual void CreateDCIOmanager(std::string, std::string) {}
      // virtual method for creating DC I/O manager for the detector

    protected:
      G4int m_verbose;

    private:
      std::string m_name;

}; // End of class G4VDCIOentry

#endif

